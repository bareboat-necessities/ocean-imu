#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
STABILITY = REPO / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_sea3_hard_spectral_driver as DRIVER  # noqa: E402


class HardSpectralDriverTests(unittest.TestCase):
    def test_contract_validates_without_surrogate_source(self):
        d = DRIVER.build()
        self.assertEqual(DRIVER.validate(d), [])
        self.assertTrue(d["HARD_SPECTRAL_DRIVER_SET_CLOSED"])
        self.assertTrue(d["COMMON_DRIVER_PHASE_PROPAGATION_CLOSED"])
        self.assertTrue(d["JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED"])
        self.assertFalse(d["FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED"])
        self.assertFalse(d["FULL_FAMILY_CORRELATED_OUTER_ENCLOSURE_CLOSED"])
        self.assertFalse(d["DEPLOYMENT_LEFT_INCLUSION_CLOSED"])

    def test_one_driver_owns_every_sample_and_channel(self):
        d = DRIVER.build()
        same = d["same_history_spectral_factorization"]
        self.assertTrue(same["same_psi_for_all_samples"])
        self.assertTrue(same["same_psi_for_translation_and_rotation"])
        self.assertTrue(same["same_psi_for_frontend_tuner_and_shipping_geometry"])
        op = d["finite_sample_operator"]
        self.assertTrue(op["one_correlated_behavior_not_cartesian_sample_boxes"])
        self.assertTrue(op["all_cross_sample_cross_axis_terms_retained"])
        self.assertIn("K=T T*", op["gram"])
        self.assertIn("K^dagger", op["membership_condition"])

    def test_forbidden_shortcuts_remain_false(self):
        d = DRIVER.build()
        for key in (
            "finite_frequency_grid_used",
            "finite_direction_grid_used",
            "finite_seeded_harmonic_generator_used",
            "trajectory_replay_used",
            "independent_sample_boxes_used",
            "independent_axis_boxes_used",
            "independent_SEA_RAO_product_used",
            "nominal_RAO_selected",
        ):
            self.assertFalse(d[key], key)
        self.assertFalse(d["driver_space"]["power_spectrum_alone_used_as_pathwise_bound"])
        self.assertFalse(d["driver_space"]["probabilistic_good_event_used"])


if __name__ == "__main__":
    unittest.main()
