"""Check the final exported study and the numbers embedded in the website."""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from notebook_studies import read_results
from verify_final_study import binomial_upper_tail,verify_source_thresholds


class FinalStudy(unittest.TestCase):
    def test_binomial_calculation(self):
        for n,count,p in [(12,3,.2),(20,0,.4),(12,12,.9),(10,7,.01)]:
            expected=sum(math.comb(n,j)*p**j*(1-p)**(n-j) for j in range(count,n+1))
            self.assertAlmostEqual(binomial_upper_tail(count,n,p),expected,places=12)

    def test_source_formulas_agree(self):
        self.assertEqual(len(verify_source_thresholds()),7)

    def test_final_counts_and_probabilities(self):
        final=ROOT/'results/final';manifest=json.loads((final/'manifest.json').read_text())
        self.assertTrue(manifest['complete']);self.assertEqual(manifest['evaluation_trials'],577800)
        self.assertEqual(manifest['calibration_matrices'],70500)
        for study in manifest['studies']:
            folder=final/study['name']
            s=read_results(folder/'summary.csv');c=read_results(folder/'calibration.csv')
            t=read_results(folder/'trials.csv.gz');a=read_results(folder/'analytic_sum_max_check.csv')
            self.assertEqual(len(t),study['evaluation_trials'])
            self.assertEqual(set(t.condition),{'null','signal'})
            self.assertTrue(s.replicates.eq(600).all());self.assertTrue(c.calibration_reps.eq(1500).all())
            self.assertTrue(a.consistent.all())
            self.assertTrue(np.allclose(a.bonferroni_cutoff,.01/1476))
        bounds=read_results(final/'thresholds/theory_threshold_checks.csv')
        self.assertEqual(len(bounds),350);self.assertTrue(bounds.consistent.all())

    def test_website_uses_canonical_notebook(self):
        import hashlib
        soup=BeautifulSoup((ROOT/'index.html').read_text(encoding='utf-8'),'html.parser')
        rows=json.loads(soup.find('script',id='empirical-results-json').string)
        meta=json.loads(soup.find('script',id='notebook-provenance-json').string)
        self.assertEqual(meta['seed'],236)
        self.assertEqual(meta['evaluation_reps'],600)
        self.assertEqual(meta['finite_matrices'],21600)
        self.assertEqual(meta['notebook_sha256'],hashlib.sha256((ROOT/meta['notebook']).read_bytes()).hexdigest())
        exported=read_results(ROOT/'results/notebook/summary.csv')
        self.assertEqual(len(rows),238)
        self.assertEqual(len(rows),len(exported))
        self.assertEqual({r['study'] for r in rows},{'finite','growing'})
        for actual,expected in zip(rows,exported.to_dict('records')):
            self.assertAlmostEqual(actual['error'],expected['error'],places=14)
            if actual['task']=='detection' and actual['study']=='finite':
                self.assertAlmostEqual(actual['error'],(actual['false_positive']+actual['false_negative'])/2)
        self.assertEqual((ROOT/'index.html').read_bytes(),(ROOT/'site/index.html').read_bytes())


if __name__=='__main__':
    unittest.main()
