#!/usr/bin/env python3
import math, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; TOOLS=ROOT/'tools'
if str(TOOLS) not in sys.path: sys.path.insert(0,str(TOOLS))
import ou3_sea3_a21_prior_free_completion as A21

class Sea3A21MixedCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.d=A21.build(); cls.failures=A21.validate(cls.d)
    def test_certificate_closes(self):
        d=self.d; self.assertEqual([],self.failures)
        self.assertEqual('COMPLETE_SEA3_NORMAL_LIVE_WORD',d['canonical_source'])
        self.assertTrue(d['A21_prior_free_completion_closed'])
        self.assertTrue(d['full_21x21_Omega_minus_delta_P_LDLT_closed'])
        self.assertTrue(d['full_21x21_interval_LDLT_used'])
        self.assertEqual(1e-18,d['useful_gate'])
    def test_no_eta9_or_full_D_shortcut(self):
        d=self.d; self.assertEqual('ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION',d['paper_active_bias_route'])
        self.assertFalse(d['eta9_point_packet_shortcut_used']); self.assertFalse(d['full_A21_D_inverse_available']); self.assertFalse(d['full_A21_prior_free_D_inverse_identity_used'])
        self.assertTrue(d['finite_tau_detectability_does_not_imply_full_A21_information_inverse'])
    def test_early_release_full_matrix(self):
        e=self.d['early']; self.assertTrue(e['covers_any_release_before_H18_margin']); self.assertTrue(e['active_A_measurement_bound_used_for_whole_word'])
        self.assertTrue(e['full_21x21_interval_LDLT_used']); self.assertTrue(e['closed']); self.assertEqual([],e['failures'])
        self.assertGreater(e['leaves'],0); self.assertGreater(e['worst_full_21x21_LDLT_pivot_lower'],0.0)
    def test_late_release_dimension_jump(self):
        l=self.d['late']; self.assertTrue(l['H18_ba_cross_zero_at_release']); self.assertTrue(l['closed']); self.assertFalse(l['another_H_word_required'])
        self.assertGreater(l['Qba_lower'],0.0); self.assertGreater(l['ba_margin_lower'],0.0)
    def test_time_varying_memory_and_rs_are_canonical(self):
        d=self.d; self.assertTrue(d['time_varying_translation_memory_full_interval_consumed']); self.assertFalse(d['translation_entrywise_lower_diagnostic_consumed'])
        self.assertFalse(d['constant_tau_over_word_assumed']); self.assertFalse(d['source_history_graph_consumed']); self.assertTrue(d['actual_applied_SpectralMSE_R_S_retained'])
        self.assertTrue(d['event_algebra_preserves_margin']); self.assertFalse(d['source_family_replaced']); self.assertFalse(d['trajectory_replay_used']); self.assertFalse(d['independent_tau_sigma_RS_source_created']); self.assertFalse(d['P3_promoted'])

if __name__=='__main__': unittest.main()
