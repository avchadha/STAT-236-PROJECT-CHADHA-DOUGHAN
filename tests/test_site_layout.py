"""Protect the user's procedure wording while checking the new branch hierarchy."""
import json
import sys
import unittest
from pathlib import Path

from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from site_layout import procedure_content,DETECTION,RECOVERY


class ProcedureLayout(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.soup=BeautifulSoup((ROOT/'index.html').read_text(encoding='utf-8'),'html.parser')
        cls.section=cls.soup.find('section',id='methods')

    def test_original_text_and_formulas_are_unchanged(self):
        before=json.loads((ROOT/'tests/fixtures/procedure_blurbs.json').read_text(encoding='utf-8'))
        self.assertEqual(before,procedure_content(self.section))

    def test_each_goal_has_its_own_ordered_branch(self):
        for goal,expected in [('detection',DETECTION),('recovery',RECOVERY)]:
            lane=self.section.find('div',id=goal+'-procedures')
            self.assertEqual([a.h3.get_text(' ',strip=True) for a in lane.select('article.method-card')],expected)
        self.assertEqual(len(self.section.select('article.method-card')),8)

    def test_boundary_then_table_end_the_section(self):
        container=self.section.find('div',class_='container')
        blocks=container.find_all(recursive=False)
        self.assertIn('inference-flow',blocks[-3].get('class',[]))
        self.assertIn('recovery-boundary',blocks[-2].get('class',[]))
        self.assertEqual(blocks[-1].get('id'),'procedure-comparison')
        self.assertNotIn('displayed below.',blocks[-1].get_text())
        self.assertEqual(len(blocks[-1].select('tbody tr')),8)

    def test_accessible_labels_and_unique_ids(self):
        ids=[t['id'] for t in self.soup.select('[id]')]
        self.assertEqual(len(ids),len(set(ids)))
        for lane in self.section.select('.procedure-lane'):
            self.assertIsNotNone(self.soup.find(id=lane['aria-labelledby']))


if __name__=='__main__':
    unittest.main()
