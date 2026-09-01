from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_p5_entrance_search_domain as E
import ou3_p5_outer_sector_capture_certificate as C


class P4P5EntranceSearchDomainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e = E.build()

    def test_validates(self):
        self.assertEqual(E.validate(self.e), [])
        self.assertEqual(self.e["P4_P5_ENTRANCE_SEARCH_DOMAIN_CERTIFICATE"], "PASS")

    def test_p5_entrance_is_45deg_so3_not_euler_box(self):
        p5 = self.e["P5_entrance"]
        self.assertEqual(p5["attitude_representation"], "SO3_GEODESIC")
        self.assertFalse(p5["componentwise_euler_box_interpretation"])
        self.assertEqual(p5["gauged_full_attitude_angle_upper_deg"], 45.0)
        self.assertLess(p5["attitude_geometry"]["cayley_norm_upper"], 1.0)
        self.assertAlmostEqual(
            p5["attitude_geometry"]["cayley_norm_upper"],
            0.8284271247461903,
            places=14,
        )

    def test_position_entrance_is_half_Hs_per_axis(self):
        p5 = self.e["P5_entrance"]
        self.assertTrue(p5["position_truth_error_bound"])
        self.assertEqual(p5["position_component_abs_error_upper_Hs_factor"], 0.5)
        self.assertGreater(p5["position_norm_upper_Hs_factor"], 0.86)
        self.assertLess(p5["position_norm_upper_Hs_factor"], 0.87)
        self.assertGreater(
            p5["Hs_below_this_m_guarantees_smaller_position_norm_than_legacy_P1_box"],
            23.0,
        )

    def test_sea_scaled_coordinates_do_not_invent_extra_hard_bounds(self):
        scaled = self.e["sea_scaled_translation_coordinates"]
        self.assertEqual(scaled["definitions"]["position"], "delta_p / Hs")
        self.assertEqual(scaled["definitions"]["velocity"], "delta_v * Ts / Hs")
        self.assertEqual(scaled["definitions"]["integral_displacement"], "delta_S / (Hs * Ts)")
        self.assertEqual(scaled["definitions"]["latent_acceleration"], "delta_a_w * Ts^2 / Hs")
        self.assertEqual(scaled["new_hard_bounds"]["position_component_abs"], "<= 0.5")
        self.assertIsNone(scaled["new_hard_bounds"]["velocity"])
        self.assertIsNone(scaled["new_hard_bounds"]["integral_displacement"])
        self.assertIsNone(scaled["new_hard_bounds"]["latent_acceleration"])
        self.assertFalse(self.e["additional_sea_scaled_v_S_aw_hard_bounds_invented"])

    def test_p4_search_ladder_is_narrower_and_improves_geometry(self):
        rows = self.e["P4_complete_word_search"]["candidate_rows"]
        self.assertEqual([r["angle_deg"] for r in rows], [30.0, 25.0, 20.0, 15.0])
        self.assertTrue(all(r["inside_45deg_entrance"] for r in rows))
        for a, b in zip(rows, rows[1:]):
            self.assertGreater(a["cayley_norm_upper"], b["cayley_norm_upper"])
            self.assertLess(a["strong_monotonicity_factor_lower"], b["strong_monotonicity_factor_lower"])
            self.assertGreater(
                a["eta_to_residual_information_ratio_upper"],
                b["eta_to_residual_information_ratio_upper"],
            )
        self.assertTrue(self.e["P4_complete_word_search"]["search_strategy_not_theorem_assumption"])

    def test_p1_box_is_not_silently_replaced(self):
        self.assertFalse(self.e["P1_conservative_handoff_box_replaced"])
        self.assertFalse(self.e["startup_propagation_of_entrance_assumed_without_proof"])

    def test_entrance_mutations_fail_closed(self):
        d = deepcopy(self.e)
        d["P5_entrance"]["componentwise_euler_box_interpretation"] = True
        self.assertNotEqual(E.validate(d), [])

        d = deepcopy(self.e)
        d["P5_entrance"]["position_component_abs_error_upper_Hs_factor"] = 0.6
        self.assertNotEqual(E.validate(d), [])

        d = deepcopy(self.e)
        d["P1_conservative_handoff_box_replaced"] = True
        self.assertNotEqual(E.validate(d), [])


class P5OuterSectorCaptureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = C.build()

    def test_validates(self):
        self.assertEqual(C.validate(self.d), [])

    def test_capture_is_immediate_and_non_microscopic(self):
        self.assertEqual(self.d["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"], "PASS")
        self.assertEqual(self.d["N_outer_words"], 0)
        self.assertGreaterEqual(self.d["outer_sector_angle_rad"], 0.80)
        self.assertTrue(self.d["all_source_handoff_branches_enter_outer_sector"])

    def test_declared_entrance_is_distinct_from_p4_complete_word_sector(self):
        e = self.d["declared_P5_entrance"]
        self.assertTrue(self.d["P5_entrance_is_distinct_from_P4_complete_word_sector"])
        self.assertTrue(self.d["P5_entrance_set_inside_outer_geometry_sector"])
        self.assertEqual(e["gauged_full_attitude_angle_upper_deg"], 45.0)
        self.assertEqual(e["position_component_abs_error_upper_Hs_factor"], 0.5)
        self.assertFalse(e["componentwise_euler_box_interpretation"])
        self.assertFalse(self.d["P1_conservative_handoff_box_replaced"])
        self.assertFalse(self.d["startup_propagation_of_entrance_assumed_without_proof"])
        self.assertFalse(self.d["P5_CAPTURE_TO_NARROWER_P4_COMPLETE_WORD_CANDIDATE_ESTABLISHED_HERE"])

    def test_p1_gauged_handoffs_fit_inside_declared_entrance(self):
        self.assertTrue(self.d["P1_source_gauged_handoffs_inside_declared_P5_entrance"])
        b = self.d["branches"]
        self.assertTrue(b["normal_gauged"]["inside_declared_P5_entrance"])
        self.assertTrue(b["timeout_gauged"]["inside_declared_P5_entrance"])
        self.assertLessEqual(
            b["normal_gauged"]["P1_handoff_cayley_norm_upper"],
            b["normal_gauged"]["declared_P5_entrance_cayley_norm_upper"],
        )
        self.assertLessEqual(
            b["timeout_gauged"]["P1_handoff_cayley_norm_upper"],
            b["timeout_gauged"]["declared_P5_entrance_cayley_norm_upper"],
        )

    def test_all_branches_are_source_faithful(self):
        b = self.d["branches"]
        self.assertTrue(b["normal_gauged"]["inside_outer_sector"])
        self.assertTrue(b["timeout_gauged"]["inside_outer_sector"])
        self.assertTrue(b["timeout_ungauged"]["inside_outer_gravity_sector"])
        self.assertTrue(b["timeout_ungauged"]["declared_entrance_tilt_inside_outer_sector"])
        self.assertFalse(b["timeout_ungauged"]["full_heading_radius_assigned"])

    def test_ungauged_timeout_uses_upper_cosine_enclosure(self):
        b = self.d["branches"]["timeout_ungauged"]
        self.assertTrue(self.d["upper_cosine_enclosure_used_for_ungauged_inclusion"])
        self.assertEqual(b["boundary_cosine_direction_used"], "UPPER_ENCLOSURE")
        self.assertLessEqual(
            self.d["outer_sector_cosine_lower"],
            self.d["outer_sector_cosine_upper"],
        )
        self.assertGreaterEqual(
            b["tilt_cosine_lower"],
            b["outer_sector_cosine_upper"],
        )

    def test_retired_microscopic_route_is_not_capture_definition(self):
        self.assertFalse(self.d["legacy_microscopic_inner_seed_used_as_outer_capture_target"])
        self.assertFalse(self.d["legacy_uniform_transport_route_used"])
        self.assertFalse(self.d["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"])

    def test_wrong_boundary_direction_or_dead_route_fails_closed(self):
        d = deepcopy(self.d)
        d["upper_cosine_enclosure_used_for_ungauged_inclusion"] = False
        self.assertNotEqual(C.validate(d), [])

        d = deepcopy(self.d)
        d["branches"]["timeout_ungauged"]["boundary_cosine_direction_used"] = "LOWER_ENCLOSURE"
        self.assertNotEqual(C.validate(d), [])

        d = deepcopy(self.d)
        d["legacy_uniform_transport_route_used"] = True
        self.assertNotEqual(C.validate(d), [])

        d = deepcopy(self.d)
        d["legacy_microscopic_inner_seed_used_as_outer_capture_target"] = True
        self.assertNotEqual(C.validate(d), [])

        d = deepcopy(self.d)
        d["P1_conservative_handoff_box_replaced"] = True
        self.assertNotEqual(C.validate(d), [])


if __name__ == "__main__":
    unittest.main()
