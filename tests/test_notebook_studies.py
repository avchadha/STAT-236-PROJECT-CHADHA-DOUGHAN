"""Coverage, data integrity, and independent mathematical checks for notebook studies."""
import ast
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments'))
from notebook_studies import (broad_plan, threshold_plan, power_law_plan, k_from_gamma,
                              reference_scales, read_results, run_study, theory_errors,
                              validate_study)


class NotebookStudies(unittest.TestCase):
    def test_original_broad_design_preserved(self):
        plan, mapping = broad_plan()
        archived = json.loads((ROOT/'tests/fixtures/broad_design.json').read_text())
        self.assertEqual(len(mapping), 35)
        self.assertEqual(len(plan), 31)
        self.assertEqual({(c['n'],c['k']) for c in plan}, {(c['n'],c['k']) for c in archived})
        for cfg in plan:
            mus = next(c['mus'] for c in archived if (c['n'],c['k'])==(cfg['n'],cfg['k']))
            np.testing.assert_allclose(sorted(cfg['mus']), sorted(mus), rtol=0, atol=1e-12)

    def test_integer_crossover(self):
        for n, gamma, expected in [(64,2/3,16), (64,.5,8), (8,2/3,4), (16,.5,4)]:
            self.assertEqual(k_from_gamma(n,gamma), expected)

    def test_all_reference_levels_bracketed(self):
        required = {'sum','scan','maximum','exact_recovery','spectral','soft_spectral'}
        measured = set()
        for cfg in threshold_plan():
            for name, scale in reference_scales(cfg['n'],cfg['k']).items():
                measured.add(name)
                signals = [s for s in cfg['signals'] if s['reference'] == name]
                self.assertEqual({s['parameter'] for s in signals}, {.5,1.,2.})
                self.assertLess(min(s['mu'] for s in signals),scale)
                self.assertGreater(max(s['mu'] for s in signals),scale)
        self.assertEqual(measured, required)
        # Independently known values, including the local-swap recovery term.
        s = reference_scales(64,4)
        self.assertEqual(s['sum'],4)
        self.assertAlmostEqual(s['scan'],math.sqrt(math.log(16)))
        self.assertAlmostEqual(s['exact_recovery'],
                               (math.sqrt(2*math.log(4))+math.sqrt(2*math.log(60)))/2)

    def test_power_law_boundaries_crossed(self):
        plan = power_law_plan()
        self.assertEqual(len(plan),12)
        for cfg in plan:
            g = math.log(cfg['k'])/math.log(cfg['n'])
            betas = [s['parameter'] for s in cfg['signals'] if s['reference'] == 'beta']
            for b in [max(g/2,2*g-1),max(0,2*g-1)]:
                self.assertTrue(any(abs(x-b) < 1e-10 for x in betas))
                self.assertTrue(any(x > b+1e-3 for x in betas))
                if b > 0:
                    self.assertTrue(any(x < b-1e-3 for x in betas))

    def test_null_csv_label_survives(self):
        data = read_results(io.StringIO('condition,value\nnull,1\nsignal,\n'))
        self.assertEqual(data.condition.tolist(),['null','signal'])
        self.assertTrue(pd.isna(data.value.iloc[1]))

    def test_short_notebook_algorithms(self):
        from functools import cache
        from itertools import combinations
        from run_experiments import recovery_spectral

        notebook = nbformat.read(ROOT/'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb',as_version=4)
        namespace = dict(np=np, cache=cache, combinations=combinations,
                         comb=math.comb, ceil=math.ceil)
        for cell in notebook.cells:
            if cell.cell_type == 'code':
                definitions = [n for n in ast.parse(cell.source).body
                               if isinstance(n, ast.FunctionDef)]
                exec(compile(ast.Module(body=definitions, type_ignores=[]),
                             '<notebook>', 'exec'), namespace)
        rng = np.random.default_rng(236)
        for n, k in [(5, 2), (6, 3)]:
            for _ in range(5):
                X = rng.normal(size=(n, n))
                value, rows, cols = namespace['scan'](X, k)
                brute = max(X[np.ix_(a, b)].sum() / k
                            for a in combinations(range(n), k)
                            for b in combinations(range(n), k))
                self.assertAlmostEqual(value, brute)
                self.assertAlmostEqual(X[np.ix_(rows, cols)].sum() / k, brute)
        for n, k, threshold in [(12, 2, 1.2), (16, 8, 0.0)]:
            X, rows, cols = namespace['make_matrix'](n, k, 3, 1, rng)
            for soft in [False, True]:
                expected = recovery_spectral(X, k, threshold=threshold if soft else None)
                actual = namespace['spectral'](X, k, threshold if soft else 0)
                for a, b in zip(actual, expected):
                    np.testing.assert_array_equal(a, b)
        # The generator changes exactly the planted entries, and nothing else.
        X, rows, cols = namespace['make_matrix'](12, 2, 3, 1, np.random.default_rng(236))
        Z = np.random.default_rng(236).standard_normal((12, 12))
        Z[np.ix_(rows, cols)] += 3
        np.testing.assert_array_equal(X, Z)
        X, rows, cols = namespace['make_matrix'](12, 2, 0, 0, np.random.default_rng(236))
        self.assertEqual(len(rows) + len(cols), 0)
        np.testing.assert_array_equal(X, np.random.default_rng(236).standard_normal((12,12)))
        self.assertEqual(namespace['scan_cutoff'](list(range(99)), .05), 94)
        self.assertEqual(namespace['recovery_error'](([1], [2]), [1], [2]), 0)
        self.assertEqual(namespace['recovery_error'](([1], [3]), [1], [2]), 1)

    def test_checkpoint_and_missing_result_detection(self):
        plan, _ = broad_plan()
        cfg = plan[0]
        with tempfile.TemporaryDirectory() as folder:
            a,t,c,_ = run_study('test',[cfg],folder,cal_reps=100,eval_reps=3)
            b,u,d,_ = run_study('test',[cfg],folder,cal_reps=100,eval_reps=3)
            pd.testing.assert_frame_equal(a,b,check_dtype=False,atol=1e-12)
            self.assertEqual(set(u.condition),{'null','signal'})
            with self.assertRaises(AssertionError):
                validate_study([cfg],a.iloc[1:],c,t,100,3)
            with self.assertRaises(AssertionError):
                validate_study([cfg],a,c,t.iloc[1:],100,3)
            bad = a.copy(); bad.loc[0,'count'] += 1
            with self.assertRaises(AssertionError):
                validate_study([cfg],bad,c,t,100,3)

    def test_theoretical_sum_error(self):
        rows = [dict(task='detection',n=16,k=4,mu=2.,condition=condition,
                     stat_sum=value, stat_maximum=0.)
                for condition in ['null','signal'] for value in [0.,3.]]
        risk = theory_errors(pd.DataFrame(rows))
        total = risk[risk.method == 'sum'].iloc[0]
        self.assertEqual(total.cutoff,1.)
        self.assertEqual(total.total_error,1.)
        self.assertAlmostEqual(total.exact_sum_risk,.31731050786291415)


if __name__ == '__main__':
    unittest.main()
