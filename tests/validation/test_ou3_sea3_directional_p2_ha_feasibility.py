from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_directional_p2_ha_feasibility as sea3  # noqa: E402


class Sea3DirectionalP2HaFeasibilityTest(unittest.TestCase):
    def test_directional_response_enclosure_retains_cross_axis_structure(self) -> None:
        d = sea3.directional_response_enclosure(ROOT)
        self.assertTrue(d["response_frequency_band_source_parity"])
        self.assertEqual(d["sea_modes_max"], 3)
        self.assertEqual(d["response_hypothesis_status"], "PROVISIONAL_NOT_DATA_DERIVED")
        self.assertFalse(d["physical_SEA0_promoted"])
        self.assertFalse(d["measured_hull_RAO_used"])
        coupling = d["directional_cross_axis_coupling"]
        self.assertTrue(coupling["rank_one_complex_outer_product_retained_before_outer_bound"])
        self.assertTrue(coupling["response_spectral_matrix_PSD"])
        self.assertFalse(coupling["independent_cartesian_axis_boxes_used"])
        self.assertFalse(d["numerical_response_gain_is_allowed_to_prune_P2"])
        self.assertGreater(
            d["matrix_moment_trace_upper_per_Hs2"]["acceleration"][1], 0.0
        )

    def test_shipping_bridge_contains_exact_clamp_and_stage_semantics(self) -> None:
        d = sea3.source_bridge_contract(ROOT)
        self.assertTrue(d["all_shipping_bridge_markers_present"])
        self.assertTrue(all(d["wrapper_markers"].values()))
        self.assertTrue(all(d["tuner_markers"].values()))

    def test_full_p2_pass_is_inherited_but_failure_is_not(self) -> None:
        inclusion = {"SEA3_TO_P2_INCLUSION_CERTIFICATE": "PASS"}
        ha = {
            "modes": {
                "H": {"relative_Riccati_injection_margin_lower": 2.0e-18},
                "A": {"relative_Riccati_injection_margin_lower": 3.0e-18},
            }
        }
        canonical_pass = {
            "P3_CANONICAL_PASS": True,
            "worst_H_A_margin": 2.0e-18,
            "useful_gate": 1.0e-18,
        }
        passed = sea3._ha_inheritance(inclusion, ha, canonical_pass)
        self.assertEqual(passed["SEA3_HA_FEASIBILITY"], "PASS_BY_P2_SUPERSET")
        self.assertTrue(passed["SEA3_HA_feasible_by_existing_uniform_certificate"])

        canonical_fail = dict(canonical_pass)
        canonical_fail["P3_CANONICAL_PASS"] = False
        failed = sea3._ha_inheritance(inclusion, ha, canonical_fail)
        self.assertEqual(
            failed["SEA3_HA_FEASIBILITY"],
            "INCONCLUSIVE_REQUIRES_SEA3_NARROWING",
        )
        self.assertFalse(failed["SEA3_HA_feasible_by_existing_uniform_certificate"])
        self.assertFalse(failed["uniform_P2_FAIL_implies_SEA3_FAIL"])


if __name__ == "__main__":
    unittest.main()
