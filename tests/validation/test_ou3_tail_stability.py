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



class RefinementTests(unittest.TestCase):
    def deployed_like(self, **kw):
        d=dict(min_samples=128,min_window_s=30.0,service_window_s=1.0,true_field_norm_lower=20.0,
               true_field_norm_upper=75.0,measurement_residual_norm=2.0,
               true_horizontal_lower=15.0,tilt_error_upper_rad=math.radians(2.0),
               max_norm_ratio_from_mean=.35,min_horizontal_fraction=.05)
        d.update(kw); return RefinementPremises(**d)

    def test_physical_residual_closes_norm_gate(self):
        p=self.deployed_like()
        m=refinement_gate_margins(p)
        self.assertLess(m["norm_ratio_upper"],.35)
        self.assertGreater(m["horizontal_fraction_lower"],.05)
        self.assertTrue(refinement_sample_gate_uniform(p))

    def test_service_closes_count_and_window_gates(self):
        p=self.deployed_like()
        # Conservative theorem bound: one usable event per one-second service
        # window. Real shipping cadence is much faster but is not needed.
        self.assertEqual(refinement_completion_bound(90.0,p),218.0)

    def test_large_tilt_error_refuses_horizontal_gate(self):
        p=self.deployed_like(tilt_error_upper_rad=math.radians(15.0))
        self.assertFalse(refinement_sample_gate_uniform(p))
        with self.assertRaises(ValueError): refinement_completion_bound(90.0,p)

    def test_residual_too_large_refuses_norm_gate(self):
        p=self.deployed_like(measurement_residual_norm=4.0)
        self.assertFalse(refinement_sample_gate_uniform(p))

if __name__=="__main__": unittest.main()
