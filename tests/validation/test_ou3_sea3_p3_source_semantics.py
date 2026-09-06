from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_p3_full_preconditions as pre  # noqa: E402
import ou3_sea3_riccati_metric_p3 as gate  # noqa: E402


class Sea3P3SourceSemanticsTest(unittest.TestCase):
    def test_preconditions_require_same_compact_phase_continuous_source(self):
        d = pre.build()
        self.assertEqual(pre.validate(d), [])
        m = d["mandatory_preconditions"]
        for key in (
            "compact_SEA3_parameter_domain_consumed",
            "compact_SEA3_transition_relation_consumed",
            "phase_continuous_SEA3_realization_required",
            "same_xs_lambda_drives_entire_source_word",
            "hard_pathwise_SEA3_conditions_retained",
            "stochastic_event_not_source_generator",
            "stochastic_event_not_homogeneous_pruner",
            "actual_applied_per_axis_R_S_required",
        ):
            self.assertTrue(m[key], key)
        self.assertFalse(d["complete_SEA3_source_family_materialized"])

    def test_canonical_gate_rejects_all_non_sea3_source_shortcuts(self):
        d = gate.build()
        self.assertEqual(gate.validate(d), [])
        self.assertTrue(d["complete_SEA3_compact_parameter_domain_consumed"])
        self.assertTrue(d["complete_SEA3_compact_transition_relation_consumed"])
        self.assertTrue(d["complete_SEA3_phase_continuous_realization_required"])
        self.assertTrue(d["same_xs_lambda_drives_entire_word"])
        self.assertTrue(d["stochastic_forcing_does_not_generate_source_words"])
        self.assertTrue(d["stochastic_forcing_does_not_prune_homogeneous_family"])
        for key in (
            "gaussian_good_event_source_used",
            "spectral_moment_only_source_used",
            "arbitrary_bounded_input_source_used",
            "independent_tau_sigma_RS_TS_extrema_product_used",
            "independent_sea_x_RAO_product_used",
            "point_source_word_used",
            "selected_four_S_word_used",
        ):
            self.assertFalse(d[key], key)
        self.assertFalse(d["P3_CANONICAL_PASS"])
        self.assertFalse(d["P4_MAY_CONSUME_P3"])


if __name__ == "__main__":
    unittest.main()
