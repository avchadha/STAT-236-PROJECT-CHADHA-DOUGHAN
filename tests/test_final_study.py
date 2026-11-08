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

    def test_website_has_only_current_results(self):
        root=ROOT/'index.html';soup=BeautifulSoup(root.read_text(encoding='utf-8'),'html.parser')
        data=json.loads(soup.find('script',id='empirical-results-json').string)
        meta=json.loads(soup.find('script',id='empirical-metadata-json').string)
        self.assertEqual({r['study'] for r in data},{'broad','thresholds','power_laws'})
        self.assertTrue(all(r['reps']==600 for r in data))
        self.assertFalse(any(r['s'] in ['D1','D2','D3','D4','R1','R2','R3','R4'] for r in data))
        self.assertEqual(meta['manifest']['evaluation_trials'],577800)
        index={(r['study'],r['task'],r['n'],r['k'],r['mu'],r['m'],r['metric']):r for r in data}
        for study in meta['manifest']['studies']:
            for r in read_results(ROOT/'results/final'/study['name']/'summary.csv').itertuples():
                key=(study['name'],r.task,r.n,r.k,r.mu,r.method,r.metric)
                self.assertAlmostEqual(index[key]['y'],r.estimate,places=14)
        self.assertEqual(root.read_bytes(),(ROOT/'site/index.html').read_bytes())


if __name__=='__main__':
    unittest.main()
