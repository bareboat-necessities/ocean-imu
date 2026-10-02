"""Synthetic algebra/negative tests; no fixture values qualify actual hardware."""
from fractions import Fraction as F
import math
import unittest
import numpy as np
from tools.stability.ou3_theorem.imu_temporal import (
    FastWindow, OpenTemporalQualification, slow_change, two_epoch_envelope,
    averaged_fast_difference, fast_weighted_outer, slow_weighted_outer,
    audit_fast_prefix,
)
from tools.stability.ou3_theorem.imu_bias import BiasLimits, BiasSample, audit_bias_trace
from tools.stability.ou3_theorem.imu_two_timescale_certificate import certificate
from tools.stability.ou3_theorem.regimes import certificate as regime_certificate
from tools.stability.ou3_theorem.finite_residual_obstruction import certificate as legacy_certificate


class ImuTemporalTests(unittest.TestCase):
    def test_slow_reachable_change_saturates_at_two_amplitudes(self):
        self.assertEqual(slow_change(.25,.001,1),.001)
        self.assertEqual(slow_change(.25,.001,1000),.5)
        self.assertEqual(slow_change(.25,.001,0),0)

    def test_raw_opposite_peaks_remain_possible_despite_fast_accumulation(self):
        w=FastWindow(1,.1)
        out=two_epoch_envelope(.25,.001,1,.01,fast=w,first_cell_s=.01,second_cell_s=.01)
        self.assertEqual(out['fast_change'],2)
        self.assertAlmostEqual(out['total_change'],2.00001)
        trace=audit_fast_prefix([0,.01,.02],[[1,0,0],[-1,0,0]],1,w,tolerance=1e-15)
        self.assertTrue(trace['finite_prefix_pass'])
        self.assertFalse(trace['all_time_certified'])

    def test_complete_cell_cap_and_identical_epoch(self):
        w=FastWindow(1,.01)
        out=two_epoch_envelope(.25,.001,1,.1,fast=w,first_cell_s=.1,second_cell_s=.2)
        self.assertAlmostEqual(out['fast_change'],.15)
        self.assertEqual(two_epoch_envelope(.25,.001,1,0,fast=w,
            first_cell_s=.1,second_cell_s=.1)['total_change'],0)
        self.assertAlmostEqual(w.cell_amplitude(2,1),.01)
        with self.assertRaises(ValueError):
            two_epoch_envelope(.25,.001,1,.01,fast=w,first_cell_s=.1,second_cell_s=.1)

    def test_unknown_temporal_profile_is_not_raw_reachability(self):
        out=two_epoch_envelope(.25,.001,1,1,fast=None,first_cell_s=.005,second_cell_s=.005)
        self.assertEqual(out['amplitude_only_outer'],2.001)
        self.assertIsNone(out['fast_change']);self.assertIsNone(out['total_change'])
        for fn in [lambda: averaged_fast_difference(None,1,1,1),
                   lambda: fast_weighted_outer([np.eye(3)],[1],1,None)]:
            with self.assertRaises(OpenTemporalQualification):fn()

    def test_window_mean_uses_overlap_and_tiled_history(self):
        w=FastWindow(2,.1)
        self.assertAlmostEqual(w.accumulation(4.05,1),.25)
        self.assertAlmostEqual(averaged_fast_difference(w,1,2,.01),.01)
        self.assertEqual(averaged_fast_difference(w,1,2,0),0)
        self.assertAlmostEqual(averaged_fast_difference(w,1,2,2),.1)

    def test_zero_final_mean_does_not_certify_all_placed_windows(self):
        # Integral over the complete record is zero; the positive half fails.
        out=audit_fast_prefix([0,.1,.2,.3,.4],[[1,0,0],[1,0,0],[-1,0,0],[-1,0,0]],
                              1,FastWindow(.4,.15))
        self.assertTrue(out['amplitude_pass']);self.assertFalse(out['finite_prefix_pass'])
        self.assertAlmostEqual(out['max_short_window_integral'],.2)

    def test_cross_word_and_regime_windows_cannot_reset(self):
        w=FastWindow(1,.15)
        for start in (0,.5):
            self.assertTrue(audit_fast_prefix([start,start+.5],[[.2,0,0]],1,w)['finite_prefix_pass'])
        joined=audit_fast_prefix([0,.5,1],[[.2,0,0],[.2,0,0]],1,w)
        self.assertFalse(joined['finite_prefix_pass'])
        lim=BiasLimits(.25,.001,.02,.00001,1,1,w,w)
        # Rate bounds alone pass; combined fast history fails in both sensors.
        samples=[BiasSample('same',t,(0,0,0),(0,0,0),(.2,0,0),(.2,0,0)) for t in (0,.5,1)]
        out=audit_bias_trace(samples,lim)
        self.assertFalse(out['finite_prefix_pass'])
        self.assertFalse(out['fast_prefixes']['accel']['finite_prefix_pass'])
        self.assertFalse(out['fast_prefixes']['gyro']['finite_prefix_pass'])

    def test_nonuniform_grid_checks_horizon_intersections_not_only_knots(self):
        t=[0,.2,.5,.9];v=[[1,0,0],[1,0,0],[-1,0,0]]
        out=audit_fast_prefix(t,v,1,FastWindow(.4,.35))
        self.assertFalse(out['finite_prefix_pass'])
        self.assertAlmostEqual(out['max_short_window_integral'],.4)

    def test_signed_matrix_weights_do_not_inherit_unweighted_cancellation(self):
        w=FastWindow(10,.2)
        vals=np.array([[.2,0,0],[-.2,0,0]])
        self.assertTrue(audit_fast_prefix([0,1,2],vals,1,w)['finite_prefix_pass'])
        self.assertAlmostEqual(np.linalg.norm(vals.sum(axis=0)),0)
        coeff=np.array([np.eye(3),-np.eye(3)])
        actual=np.linalg.norm(np.einsum('kij,kj->i',coeff,vals))
        self.assertAlmostEqual(actual,.4)
        self.assertGreaterEqual(fast_weighted_outer(coeff,[1,1],1,w),actual-1e-15)
        # A constant coefficient still has an endpoint term.
        self.assertAlmostEqual(fast_weighted_outer([np.eye(3)],[1],1,w),.2)

    def test_signed_slow_weights_retain_predecessor_rate(self):
        self.assertAlmostEqual(slow_weighted_outer([-np.eye(3),np.eye(3)],[1],.25,.001),.001)
        self.assertAlmostEqual(slow_weighted_outer([-np.eye(3),np.eye(3)],[1000],.25,.001),.5)

    def test_finite_prefix_never_promotes_to_all_time(self):
        out=audit_fast_prefix([0,1],[[0,0,0]],1,FastWindow(2,.1))
        self.assertTrue(out['finite_prefix_pass']);self.assertFalse(out['all_time_certified'])
        self.assertIsNone(audit_fast_prefix([0,1],[[0,0,0]],1,None)['finite_prefix_pass'])

    def test_both_witnesses_have_the_correct_new_model_status(self):
        report=certificate(); w=report['oscillatory_witness']
        self.assertEqual(F(w['accel_slow_required_rate_mps3']),F('.04903325'))
        self.assertEqual(F(w['gyro_slow_required_rate_rad_s2']),F('.0025'))
        self.assertEqual(F(w['gyro_rate_to_candidate_limit_ratio']),250)
        self.assertEqual(w['current_model_admissibility'],'OPEN')
        self.assertFalse(w['all_slow_admissible_under_candidate_rates'])
        self.assertGreater(F(w['delivered_sample_necessary_accel_K_lower_mps']),F('.37065'))
        self.assertGreater(F(w['delivered_sample_necessary_gyro_K_lower_rad']),F('.019708'))
        r=regime_certificate()
        self.assertTrue(r['fast_accel_and_gyro_identically_zero'])
        self.assertTrue(all(F(v)>0 for v in r['witness_bias_margin'].values()))
        self.assertGreater(F(r['witness_actual_service_lower']),1)
        legacy=legacy_certificate()
        self.assertTrue(legacy['finite_residual_inner_retention_refuted_on_this_class'])
        self.assertIn('OPEN',legacy['current_two_timescale_admissibility'])
        self.assertFalse(legacy['current_model_storage_lower_bound_inferred'])
        self.assertFalse(report['marine_and_magnetic']['stronger_EXCITED_MOVING_added'])

    def test_invalid_profiles_fail_closed(self):
        for value in (-1,math.inf,math.nan):
            with self.assertRaises(ValueError):FastWindow(value,.1)
            with self.assertRaises(ValueError):FastWindow(1,value)
            with self.assertRaises(ValueError):slow_change(1,.1,value)
        with self.assertRaises(ValueError):FastWindow(0,0)
        with self.assertRaises(ValueError):FastWindow(1,1).validate_amplitude(1)
        with self.assertRaises(ValueError):FastWindow(1,1).validate_amplitude(0)
        self.assertEqual(FastWindow(1,0).accumulation(5,0),0)
        with self.assertRaises(ValueError):audit_fast_prefix([0,0],[[0,0,0]],1,None)
        with self.assertRaises(ValueError):fast_weighted_outer([np.eye(3)],[0],1,FastWindow(1,.1))


if __name__=='__main__':unittest.main()

class LinkedImuSupplyTests(unittest.TestCase):
    def test_both_profiles_required_and_matrix_action_retained(self):
        import numpy as np
        from tools.stability.ou3_theorem.imu_bias import BiasLimits
        from tools.stability.ou3_theorem.linked_supply import imu_supply_outer
        ops={k:np.array([np.eye(3),-np.eye(3)])*.5 for k in
             ('accel_slow','accel_fast','gyro_slow','gyro_fast')}
        from dataclasses import replace
        open_limits=replace(BiasLimits.from_constants(),accel_fast_window=None,gyro_fast_window=None)
        with self.assertRaises(OpenTemporalQualification):
            imu_supply_outer(ops,[.5,.5],open_limits)
        limits=BiasLimits(1,.01,1,.01,1,1,FastWindow(1,.1),FastWindow(1,.1))
        out=imu_supply_outer(ops,[.5,.5],limits)
        self.assertAlmostEqual(out['terms']['accel_slow'],.0025)
        self.assertAlmostEqual(out['terms']['accel_fast'],.2)
        self.assertFalse(out['linked_chi_supremum_verified'])
        self.assertFalse(out['physical_device_qualified'])
        with self.assertRaises(ValueError):
            imu_supply_outer({'accel_fast':ops['accel_fast']},[.5,.5],limits)

if __name__ == "__main__":
    unittest.main()
