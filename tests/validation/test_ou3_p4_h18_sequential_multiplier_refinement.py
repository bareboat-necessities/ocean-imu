#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_sequential_multiplier_refinement as R

class TestH18SequentialMultiplierRefinement(unittest.TestCase):
    def test_refinement_executes_fail_closed(self):
        d=R.build();self.assertEqual(R.validate(d),[])
        self.assertTrue(d['three_physical_sector_multipliers_refined_independently'])
        self.assertIn('candidate_found',d['refinement'])
        self.assertFalse(d['production_endpoint_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])
    def test_neighborhood_is_positive_and_asymmetric_capable(self):
        v=R._near(2.0);self.assertTrue(all(x>0 for x in v));self.assertGreater(len(v),3)
        self.assertNotEqual(v[0],v[-1])
if __name__=='__main__':unittest.main()
