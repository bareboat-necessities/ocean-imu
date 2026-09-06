#!/usr/bin/env python3
import math, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; TOOLS=ROOT/'tools'
if str(TOOLS) not in sys.path: sys.path.insert(0,str(TOOLS))
import ou3_sea3_a21_prior_free_completion as A21

class Sea3A21PriorFreeCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.d=A21.build(); cls.failures=A21.validate(cls.d)
    def test_certificate_closes(self):
        d=self.d; self.assertEqual([],self.failures)
        self.assertEqual('COMPLETE_SEA3_NORMAL_LIVE_WORD',d['canonical_source'])
        self.assertTrue(d['A21_prior_free_completion_closed'])
        self.assertTrue(d['full_21x21_Omega_minus_delta_P_LDLT_closed'])
        self.assertTrue(d['full_21x21_interval_LDLT_used'])
        self.assertEqual(1e-18,d['useful_gate'])
        self.assertGreater(d['x_cells_certified'],0)
        self.assertEqual([],d['x_cell_failures'])
        self.assertGreater(d['worst_full_A21_LDLT_pivot_lower'],0.0)
    def test_no_eta9_and_full_information_is_proved_by_finite_estimator(self):
        d=self.d
        self.assertEqual('ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION',d['paper_active_bias_route'])
        self.assertFalse(d['eta9_point_packet_shortcut_used'])
        self.assertTrue(d['finite_full_A21_linear_estimator_constructed'])
        self.assertTrue(d['finite_full_A21_estimator_implies_D_strictly_positive'])
        self.assertTrue(d['full_A21_prior_free_D_inverse_identity_used'])
    def test_A_mode_endpoint_accelerometer_completes_ba(self):
        d=self.d
        self.assertTrue(d['H18_finite_memory_estimator_consumed'])
        self.assertTrue(d['H18_estimator_pays_active_ba_nuisance'])
        self.assertTrue(d['required_A_mode_accelerometer_row_consumed'])
        self.assertTrue(d['A_mode_accelerometer_H_ba_is_identity'])
        e=d['A21_diffuse_prior_covariance_upper']
        self.assertGreater(e['H18_trace_upper'],0.0)
        self.assertGreater(e['ba_estimator_trace_upper'],0.0)
        self.assertGreater(e['full_A21_trace_upper'],e['H18_trace_upper'])
        self.assertTrue(e['no_H_estimator_measurement_noise_independence_assumed'])
    def test_shipping_process_and_event_algebra_are_retained(self):
        d=self.d
        self.assertTrue(d['stable_factored_shipping_integrated_OU_Q_consumed'])
        self.assertTrue(d['shipping_active_ba_GM_Q_consumed'])
        self.assertTrue(d['actual_applied_SpectralMSE_R_S_retained_through_H18_component'])
        self.assertTrue(d['event_algebra_preserves_margin_after_closure'])
        self.assertFalse(d['old_one_step_Euclidean_Q_min_used'])
        self.assertFalse(d['scalar_beta_contraction_used'])
        self.assertFalse(d['blockwise_minimum_contraction_used'])
    def test_no_replacement_source(self):
        d=self.d
        for k in ('source_family_replaced','trajectory_replay_used','independent_tau_sigma_RS_source_created','P3_promoted'):
            self.assertFalse(d[k],k)

if __name__=='__main__': unittest.main()
