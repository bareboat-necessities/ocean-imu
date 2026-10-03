from fractions import Fraction
from dataclasses import replace
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.imu_temporal import FastWindow
from tools.stability.ou3_theorem.imu_bias import (
    BiasContinuationCertificate,BiasLimits,BiasSample,accel_prediction_error,
    accel_prediction_relation,audit_bias_trace,continuation_admitted,
    correction_error,estimator_prediction_factor,projected_error,
    projection_defect,release_error,successor_allowed,project_estimate,
    prediction_joint_blocks,correction_projection_relation,
)

LIMITS=BiasLimits(0.22516660498395405,0.001,0.02,1.0e-5,.3,.02,
                  FastWindow(2,.1),FastWindow(2,.01))

class ImuBiasTests(unittest.TestCase):
    def test_successor_is_not_an_independent_box(self):
        self.assertFalse(successor_allowed((0,0,0),(0.1,0,0),0.005,LIMITS.B_a_s_mps2,LIMITS.D_a_s_mps3))

    def test_same_history_slow_drift_passes(self):
        s=[BiasSample("h",0.0,(0.01,0,0),(0.001,0,0)),
           BiasSample("h",1.0,(0.0105,0,0),(0.001005,0,0))]
        self.assertTrue(audit_bias_trace(s,LIMITS)["finite_prefix_pass"])

    def test_all_time_slow_fast_certificate_is_fail_closed(self):
        # Synthetic conditional certificate, not an assembled-device claim.
        good=BiasContinuationCertificate(
            "h",.20,5e-4,.01,5e-6,True,True,True,.2,.01,
            FastWindow(2,.1),FastWindow(2,.01),"split-1","synthetic-proof",
            True,True)
        self.assertTrue(continuation_admitted(good,LIMITS))
        for kwargs in (
            dict(predecessor_continuity_certified=False),
            dict(calibration_scope_qualified=False),
            dict(all_time_continuation_certified=False),
            dict(accel_fast_window=None),dict(gyro_fast_window=None),
            dict(all_placed_windows_certified=False),
            dict(delivered_sample_scope_certified=False),
            dict(decomposition_id=""),dict(temporal_evidence_id=""),
            dict(accel_fast_window=FastWindow(2,.101)),
            dict(gyro_fast_window=FastWindow(3,.01)),
            dict(accel_fast_norm_upper_mps2=.31),
        ):
            with self.subTest(kwargs=kwargs):
                self.assertFalse(continuation_admitted(replace(good,**kwargs),LIMITS))
        self.assertFalse(continuation_admitted(good,BiasLimits.from_constants()))

    def test_missing_fast_parameters_do_not_mean_independent_noise(self):
        limits=replace(BiasLimits.from_constants(),accel_fast_window=None,gyro_fast_window=None)
        self.assertFalse(limits.temporal_parameters_present)
        out=audit_bias_trace([BiasSample("h",0,(0,0,0),(0,0,0)),
                              BiasSample("h",1,(0,0,0),(0,0,0))],limits)
        self.assertIsNone(out["finite_prefix_pass"])
        self.assertTrue(out["slow_and_amplitude_prefix_pass"])
        self.assertEqual(out["temporal_qualification"],"OPEN")
        self.assertFalse(out["all_time_certified"])

    def test_cannot_reselect_split_at_release(self):
        samples=[BiasSample("h",0,(0,0,0),(0,0,0),decomposition_id="H18"),
                 BiasSample("h",1,(0,0,0),(0,0,0),decomposition_id="A21")]
        self.assertFalse(audit_bias_trace(samples,LIMITS)["finite_prefix_pass"])

    def test_one_first_class_relation_covers_held_and_active(self):
        phi=math.exp(-0.1/5000.0)
        self.assertEqual(estimator_prediction_factor("H18",phi),1.0)
        self.assertEqual(estimator_prediction_factor("A21",phi),phi)
        held=accel_prediction_relation("H18",phi,(0.02,0,0),(0.1,0,0),(0.001,0,0))
        active=accel_prediction_relation("A21",phi,(0.02,0,0),(0.1,0,0),(0.001,0,0))
        self.assertEqual(held.phi_e,1.0)
        self.assertEqual(held.model_mismatch,(0.0,0.0,0.0))
        self.assertAlmostEqual(held.e_b_minus_next[0],0.021)
        self.assertEqual(active.phi_e,phi)
        self.assertAlmostEqual(active.model_mismatch[0],(1.0-phi)*0.1)
        expected=phi*0.02+(1.0-phi)*0.1+0.001
        self.assertAlmostEqual(active.e_b_minus_next[0],expected)
        self.assertEqual(
            active.e_b_minus_next,
            accel_prediction_error("A21",phi,(0.02,0,0),(0.1,0,0),(0.001,0,0)),
        )

    def test_correction_release_and_projection_are_separate_operations(self):
        e_corr=correction_error((0.05,0,0),(0.01,0,0))
        self.assertEqual(e_corr,(0.04,0.0,0.0))
        self.assertEqual(release_error(e_corr),e_corr)
        defect=projection_defect((0.5,0,0),0.4)
        self.assertAlmostEqual(defect[0],0.1)
        e_plus=projected_error((0.2,0,0),(0.5,0,0),0.4)
        self.assertAlmostEqual(e_plus[0],-0.2)

    def test_joint_prediction_has_one_shared_physical_increment(self):
        # Dyadic operands make this an exact rational check, not a tolerance fit.
        for mode in ("H18", "A21"):
            for phi in (0.0, 0.5, 0.875, 1.0):
                transition, driver = prediction_joint_blocks(mode, phi)
                self.assertEqual(driver, ((1.0,), (1.0,)))
                e, b, w = Fraction(1, 8), Fraction(1, 4), Fraction(1, 1024)
                actual_phi = Fraction(1 if mode == "H18" else phi)
                out = [sum(Fraction(c) * z for c, z in zip(row, (e, b)))
                       + Fraction(driver[i][0]) * w
                       for i, row in enumerate(transition)]
                estimate_next = actual_phi * (b - e)
                self.assertEqual(out[0], b + w - estimate_next)
                self.assertEqual(out[1], b + w)

    def test_projection_branch_matches_disabled_radius_and_finite_boundary(self):
        estimate = (0.5, -0.25, 0.125)
        for radius in (0.0, -0.1):
            self.assertEqual(project_estimate(estimate, radius), estimate)
            self.assertEqual(projection_defect(estimate, radius), (0.0, 0.0, 0.0))
        self.assertEqual(project_estimate((0.4, 0.0, 0.0), 0.4), (0.4, 0.0, 0.0))
        self.assertEqual(project_estimate((float("nan"), 0.0, 0.0), 0.4), (0.0, 0.0, 0.0))

    def test_complete_bias_event_retains_correction_and_projection_defect(self):
        rel = correction_projection_relation((0.2, 0, 0), (0.35, 0, 0), (0.15, 0, 0), 0.4)
        for i in range(3):
            self.assertAlmostEqual(rel.e_corrected[i], rel.e_minus[i] - rel.applied_increment[i])
            self.assertAlmostEqual(rel.e_plus[i], rel.e_corrected[i] + rel.projection_defect[i])
        self.assertEqual(rel.b_true, (0.2, 0.0, 0.0))
        self.assertAlmostEqual(math.hypot(*rel.estimate_plus), 0.4)
        # The Euclidean projection sector requires truth inside the ball.
        sector = sum(e * d for e, d in zip(rel.e_plus, rel.projection_defect))
        self.assertLessEqual(sector, 0.0)
        self.assertLessEqual(sum(e*e for e in rel.e_plus),
                             sum(e*e for e in rel.e_corrected) - sum(d*d for d in rel.projection_defect))
        # This lemma must not be claimed when the physical bias lies outside.
        outside = correction_projection_relation((0.8,0,0), (0.35,0,0), (0.15,0,0), 0.4)
        self.assertGreater(sum(e*d for e,d in zip(outside.e_plus,outside.projection_defect)), 0)

    def test_rate_admission_rejects_nonfinite_limits_and_tolerances(self):
        for invalid in (math.inf, math.nan, -0.01):
            self.assertFalse(successor_allowed((0,0,0),(0.1,0,0),0.005,invalid,0.001))
            self.assertFalse(successor_allowed((0,0,0),(0.1,0,0),0.005,1,invalid))
            self.assertFalse(successor_allowed((0,0,0),(0.1,0,0),0.005,1,0.001,tolerance=invalid))
            with self.assertRaises(ValueError):
                audit_bias_trace([BiasSample("h",0,(0,0,0),(0,0,0))], LIMITS, tolerance=invalid)

if __name__=="__main__":
    unittest.main()
