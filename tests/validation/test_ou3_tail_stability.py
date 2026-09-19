from __future__ import annotations
import math,sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.tail_stability import *

class TailStabilityTests(unittest.TestCase):
    def test_h18_is_finite_bridge(self):
        self.assertAlmostEqual(finite_bridge_bound(2,3,1,.5),3.5)
        self.assertAlmostEqual(finite_bridge_bound(2,2,2,1),11)
    def test_service_makes_count_guard_finite_after_reference_release(self):
        self.assertEqual(release_bound_from_service(250,1,1,12),263)
    def test_nonlinear_small_gain(self):
        p=A21TailPremises(.81,.05,2,.5,3)
        self.assertAlmostEqual(nonlinear_tail_ratio(p),.9025)
        self.assertTrue(math.isfinite(practical_radius_bound(p,.1)))
    def test_noncontractive_tail_refused(self):
        p=A21TailPremises(.99,.01,1,1,1)
        self.assertGreaterEqual(nonlinear_tail_ratio(p),1)
        with self.assertRaises(ValueError): practical_radius_bound(p,.1)
    def test_route_fail_closed(self):
        r=proof_route_status(reference_refinement_finite=True,h18_bridge_retained=True,
            a21_linear_uniform=True,a21_nonlinear_bound=False,a21_prefix_retained=True)
        self.assertEqual(r["h18_role"],"finite bridge only"); self.assertFalse(r["route_closed"])
    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): finite_bridge_bound(1,-1,1,0)
        with self.assertRaises(ValueError): release_bound_from_service(-1,1,1,1)
        with self.assertRaises(ValueError): A21TailPremises(-1,0,1,1,1)

if __name__=="__main__": unittest.main()
