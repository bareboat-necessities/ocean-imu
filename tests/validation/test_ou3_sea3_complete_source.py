from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_complete_source as mod  # noqa: E402


class Sea3CompleteSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.failures = mod.validate(cls.d)

    def test_complete_source_contract_is_valid_but_not_materialized(self):
        self.assertEqual(self.failures, [])
        self.assertTrue(self.d["P3_source_contract_ready"])
        self.assertFalse(self.d["P3_source_family_materialized"])
        self.assertFalse(self.d["P3_promoted"])

    def test_sea3_is_compact_and_phase_continuous(self):
        sea = self.d["SEA3_surface_family"]
        self.assertTrue(sea["parameter_domain_compact"])
        self.assertTrue(sea["compact_transition_relation_is_theorem_domain"])
        dyn = self.d["SEA3_dynamic_realization"]
        self.assertTrue(dyn["phase_continuous"])
        self.assertEqual(dyn["shaping_state"], "x^s")
        self.assertEqual(dyn["augmented_source_state"], "zeta=(x^s,lambda,z^t,q)")
        self.assertTrue(
            dyn["same_realization_drives_translation_rotation_frontend_tuner_geometry"]
        )
        self.assertFalse(dyn["probabilistic_event_may_substitute_for_realization"])
        self.assertFalse(dyn["arbitrary_bounded_input_may_substitute_for_realization"])

    def test_no_shortcut_source_can_promote(self):
        self.assertTrue(self.d["no_fallback_generators"])
        self.assertTrue(all(v is False for v in self.d["no_fallback_generators"].values()))
        stochastic = self.d["stochastic_forcing_corollary"]
        self.assertFalse(stochastic["used_to_generate_P3_source_words"])
        self.assertFalse(stochastic["used_to_prune_homogeneous_P3_family"])
        coupling = self.d["SEA3_response_couplings"]
        self.assertTrue(
            coupling["only_same_phase_continuous_SEA3_realization_may_generate_P3_words"]
        )
        self.assertTrue(coupling["moment_or_probability_bound_may_not_generate_P3_word"])

    def test_rs_remains_actual_regularizer_in_full_word(self):
        rs = self.d["R_S_regularizer"]
        self.assertEqual(rs["deployed_law"], "SpectralMSE")
        self.assertTrue(rs["actual_applied_R_S_required_at_every_due_S_update"])
        self.assertTrue(rs["all_due_S_updates_remain_in_full_word"])
        self.assertTrue(rs["full_P_column_S_cross_covariance_action_required"])
        self.assertTrue(rs["R_S_may_not_be_replaced_by_process_strictness"])


if __name__ == "__main__":
    unittest.main()
