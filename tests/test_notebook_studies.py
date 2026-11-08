"""Coverage, data integrity, and independent mathematical checks for notebook studies."""
import ast
import io
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
        archived = read_results(ROOT / 'results/gamma_sweep/summary.csv')
        self.assertEqual(len(mapping), 35)
        self.assertEqual(len(plan), 31)
        self.assertEqual({(c['n'],c['k']) for c in plan}, set(zip(archived.n,archived.k)))
        for cfg in plan:
            mus = archived[(archived.n == cfg['n']) & (archived.k == cfg['k'])].mu.unique()
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

    def test_notebook_embeds_identical_algorithms(self):
        notebook = nbformat.read(ROOT/'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb',as_version=4)
        functions = {}
        for cell in notebook.cells:
            if cell.cell_type == 'code':
                functions.update({n.name:ast.dump(n, include_attributes=False)
                                  for n in ast.parse(cell.source).body if isinstance(n,ast.FunctionDef)})
        source = ast.parse((ROOT/'experiments/run_experiments.py').read_text())
        for node in source.body:
            if isinstance(node,ast.FunctionDef) and node.name not in {'main','export_csv'}:
                self.assertEqual(functions[node.name],ast.dump(node,include_attributes=False),node.name)

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
