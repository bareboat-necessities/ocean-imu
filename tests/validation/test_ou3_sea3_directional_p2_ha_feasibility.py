from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_directional_p2_ha_feasibility as sea3  # noqa: E402


class Sea3DirectionalP2HaFeasibilityTest(unittest.TestCase):
    def test_directional_response_covers_a_continuum_not_one_rao(self) -> None:
        d = sea3.directional_response_enclosure(ROOT)
        self.assertTrue(d["shipping_sigma_band_source_parity"])
        self.assertEqual(d["sea_modes_max"], 3)
        self.assertEqual(d["response_hypothesis_status"], "ROBUST_CONTINUUM_ENVELOPE_FAMILY")
        self.assertFalse(d["physical_SEA0_promoted"])
        self.assertFalse(d["single_nominal_RAO_used"])
        self.assertFalse(d["finite_RAO_grid_used"])
        self.assertTrue(d["six_dof_parent_RAO_allowed"])
        self.assertTrue(d["single_worst_envelope_proves_entire_parameter_box_by_monotonicity"])

        box = d["rao_envelope_parameter_box"]
        self.assertEqual(box["peak_translation_gain"], [0.0, 4.0])
        self.assertEqual(box["rolloff_corner_hz"], [0.03, 1.2])
        self.assertEqual(box["high_frequency_rolloff_power_min"], 2.0)
        self.assertEqual(box["complex_phase"], "arbitrary")
        self.assertEqual(box["heading_dependence"], "arbitrary")

        coupling = d["directional_cross_axis_coupling"]
        self.assertTrue(coupling["rank_one_complex_outer_product_retained_before_outer_bound"])
        self.assertTrue(coupling["response_spectral_matrix_PSD"])
        self.assertTrue(coupling["arbitrary_phase_is_covered"])
        self.assertFalse(coupling["independent_cartesian_axis_boxes_used"])

    def test_worst_envelope_corner_dominates_every_checked_family_member(self) -> None:
        d = sea3.directional_response_enclosure(ROOT)
        worst = d["worst_envelope_trace_upper_per_Hs2"]
        samples = (
            (0.0, 0.03, 2.0),
            (0.5, 0.1, 8.0),
            (1.0, 0.2, 4.0),
            (2.0, 0.6, 3.0),
            (3.0, 1.0, 2.1),
            (4.0, 1.2, 2.0),
        )
        for gain, corner, power in samples:
            member = sea3.evaluate_rao_envelope_member(gain, corner, power, d)
            for key in ("displacement", "velocity", "acceleration"):
                self.assertLessEqual(member[key][1], worst[key][1])
                self.assertTrue(math.isfinite(member[key][1]))

        with self.assertRaises(ValueError):
            sea3.evaluate_rao_envelope_member(4.01, 1.2, 2.0, d)
        with self.assertRaises(ValueError):
            sea3.evaluate_rao_envelope_member(4.0, 1.21, 2.0, d)
        with self.assertRaises(ValueError):
            sea3.evaluate_rao_envelope_member(4.0, 1.2, 1.99, d)

    def test_rolloff_closes_unbanded_acceleration_moment(self) -> None:
        d = sea3.directional_response_enclosure(ROOT)
        theorem = d["uniform_moment_theorem"]
        self.assertTrue(theorem["analytical_not_sampled"])
        self.assertTrue(theorem["unbanded_acceleration_moment_finite"])
        self.assertEqual(theorem["proof_corner"], "G=G_max, f_c=f_c,max, p=2")
        # Replacing the old flat 4x gain at the 6-Hz band endpoint by the
        # certified p>=2 response envelope tightens the acceleration-moment
        # coefficient by essentially (6/1.2)^4 = 625.
        self.assertGreater(
            d["acceleration_moment_tightening_vs_flat_6Hz_corner_lower"],
            624.0,
        )

    def test_shipping_bridge_contains_exact_clamp_and_stage_semantics(self) -> None:
        d = sea3.source_bridge_contract(ROOT)
        self.assertTrue(d["all_shipping_bridge_markers_present"])
        self.assertTrue(all(d["wrapper_markers"].values()))
        self.assertTrue(all(d["tuner_markers"].values()))

    def test_full_p2_pass_is_inherited_for_the_whole_rao_box_but_failure_is_not(self) -> None:
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
        self.assertTrue(passed["uniform_over_entire_RAO_parameter_box"])
        self.assertEqual(passed["unchanged_useful_gate"], 1.0e-18)

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
