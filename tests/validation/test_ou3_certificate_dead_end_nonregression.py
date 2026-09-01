from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_source_path_reachability as PATH
import ou3_p4_operation_matched_sector_certificate as P4
import ou3_p5_outer_sector_capture_certificate as P5
import ou3_p4_post_translation_bottleneck as POST

DOMAIN = ROOT / "tools" / "ou3_proof_operating_domain.json"
WORKFLOW = ROOT / ".github" / "workflows" / "ou3-usable-certificates-fast.yml"


class OU3CertificateDeadEndNonRegressionTests(unittest.TestCase):
    """Keep already-disproved proof decompositions from returning silently.

    These are theorem non-regression tests, not style preferences.  PR #441
    proved that the scalar uniform P4->P5 transport route has a milliradian hard
    ceiling.  PR #438 supplied valuable complete-word translation machinery but
    also demonstrated that a widened translation block is only partial evidence
    until the complete full-state cross coupling and nonlinear return map close.
    """

    @classmethod
    def setUpClass(cls):
        cls.path = PATH.build()
        cls.p4 = P4.build()
        cls.p5 = P5.build()
        cls.domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
        cls.workflow = WORKFLOW.read_text(encoding="utf-8")

    def test_declared_deployment_domain_cannot_be_shrunk_to_make_certificate_pass(self):
        d = self.domain
        self.assertIs(d.get("trajectory_fit"), False)

        s = d["startup"]
        self.assertGreaterEqual(s["initial_non_gravitational_specific_force_norm_upper_mps2"], 4.0)
        self.assertGreaterEqual(s["initial_tangent_gyro_bias_norm_upper_rad_s"], 0.01)
        self.assertGreaterEqual(s["world_averaged_gravity_direction_error_upper_rad"], 0.02)
        self.assertGreaterEqual(s["mahony_chart_theta_star_deg"], 60.0)
        self.assertGreaterEqual(s["internal_heading_gauge_error_upper_rad"], 0.17453292519943295)
        handoff = s["physical_handoff_coordinate_bounds"]
        floors = {
            "gyro_bias_error_norm_upper_rad_s": 0.01,
            "velocity_error_norm_upper_mps": 5.0,
            "position_error_norm_upper_m": 20.0,
            "integral_displacement_error_norm_upper_m_s": 300.0,
            "accelerometer_bias_error_norm_upper_mps2": 0.5,
        }
        for key, baseline in floors.items():
            self.assertGreaterEqual(handoff[key], baseline, key)

        # The latent-acceleration startup envelope is a separately reviewed
        # physical theorem assumption.  It was corrected from the old near-1g
        # placeholder to 0.3g for wave startup; freeze that declared value so a
        # later proof cannot silently shrink it further just to pass a gate.
        expected_aw = 0.3 * s["gravity_mps2"]
        self.assertEqual(s["latent_acceleration_error_fraction_g"], 0.3)
        self.assertEqual(handoff["latent_acceleration_error_norm_upper_mps2"], expected_aw)
        self.assertEqual(
            s["latent_acceleration_error_bound_role"],
            "DECLARED_STARTUP_WAVE_ERROR_ENVELOPE_0P3G_NOT_REPLAY_FIT_OR_PROOF_RETUNING",
        )

        live = d["normal_live"]
        self.assertLessEqual(live["specific_force_norm_lower_mps2"], 5.0)
        self.assertLessEqual(live["magnetic_vector_norm_lower_uT"], 10.0)
        self.assertLessEqual(live["vector_sine_separation_lower"], 0.1)
        self.assertGreaterEqual(live["body_rate_norm_upper_deg_s"], 30.0)
        self.assertGreaterEqual(live["vector_pe_recurrence_window_s"], 1.0)
        self.assertGreaterEqual(live["active_accelerometer_bias_state_norm_upper_mps2"], 0.45)

        q = live["gravity_quotient"]
        self.assertGreaterEqual(q["accepted_packet_recurrence_window_s"], 1.0)
        self.assertLessEqual(q["accepted_packet_separation_min_s"], 0.04)
        self.assertGreaterEqual(q["accepted_packet_separation_max_s"], 1.0)
        self.assertGreaterEqual(q["non_gravitational_specific_force_norm_upper_mps2"], 4.0)

    def test_source_path_graph_cannot_omit_shipping_tuner_branches(self):
        self.assertEqual(PATH.validate(self.path), [])
        self.assertEqual(self.path["P2_SOURCE_PATH_CERTIFICATE"], "PASS")
        self.assertTrue(self.path["path_graph_ready"])
        self.assertGreater(self.path["transition_edges"], 0)
        self.assertGreater(self.path["old_worst_corner_state_count"], 0)
        self.assertFalse(self.path["usable_P4_promoted"])

        self.assertTrue(self.path["raw_tuner_sigma_subfloor_states_included"])
        self.assertLess(
            self.path["raw_tuner_sigma_partition_lower"],
            self.path["filter_sigma_floor_mps2"],
        )
        self.assertTrue(self.path["filter_sigma_floor_separate_from_tuner_state"])
        self.assertTrue(self.path["validated_exponential_used_for_ema"])
        self.assertTrue(self.path["arbitrary_late_commit_overapproximated"])
        self.assertTrue(self.path["RS_discrepancy_slew_horizon_covered"])
        self.assertTrue(self.path["RS_target_full_deployed_clamp_overapprox"])
        self.assertFalse(self.path["RS_target_powf_tightening_used"])

        bad = deepcopy(self.path)
        bad["usable_P4_promoted"] = True
        self.assertTrue(PATH.validate(bad))

    def test_p3_relative_riccati_delta_cannot_be_used_as_state_amplitude_radius(self):
        self.assertEqual(P4.validate(self.p4), [])
        self.assertIs(self.p4["whole_word_weakest_P3_delta_used_as_attitude_sector_margin"], False)
        bad = deepcopy(self.p4)
        bad["whole_word_weakest_P3_delta_used_as_attitude_sector_margin"] = True
        self.assertTrue(P4.validate(bad))

    def test_global_packet_count_times_scalar_lipschitz_defect_is_retired(self):
        self.assertIs(self.p4["global_packet_count_times_lipschitz_defect_used"], False)
        bad = deepcopy(self.p4)
        bad["global_packet_count_times_lipschitz_defect_used"] = True
        self.assertTrue(P4.validate(bad))

    def test_p4_may_not_regress_to_a_microscopic_sector(self):
        self.assertGreaterEqual(self.p4["design_full_attitude_angle_rad"], 0.80)
        bad = deepcopy(self.p4)
        bad["design_full_attitude_angle_rad"] = 0.00124
        self.assertTrue(P4.validate(bad))

    def test_operation_matched_sector_is_not_complete_word_p4_by_itself(self):
        self.assertIs(self.p4["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"], False)
        bad = deepcopy(self.p4)
        bad["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"] = True
        self.assertTrue(P4.validate(bad))

    def test_p5_cannot_return_to_legacy_uniform_transport_or_inner_seed_target(self):
        self.assertEqual(P5.validate(self.p5), [])
        self.assertIs(self.p5["legacy_uniform_transport_route_used"], False)
        self.assertIs(self.p5["legacy_microscopic_inner_seed_used_as_outer_capture_target"], False)

        bad = deepcopy(self.p5)
        bad["legacy_uniform_transport_route_used"] = True
        self.assertTrue(P5.validate(bad))

        bad = deepcopy(self.p5)
        bad["legacy_microscopic_inner_seed_used_as_outer_capture_target"] = True
        self.assertTrue(P5.validate(bad))

    def test_ungauged_timeout_uses_upper_cosine_enclosure_and_no_heading_radius(self):
        u = self.p5["branches"]["timeout_ungauged"]
        self.assertIs(u["full_heading_radius_assigned"], False)
        self.assertEqual(u["attitude_representation"], "GRAVITY_DIRECTION_QUOTIENT")
        self.assertEqual(u["boundary_cosine_direction_used"], "UPPER_ENCLOSURE")
        self.assertGreaterEqual(u["tilt_cosine_lower"], u["outer_sector_cosine_upper"])

        bad = deepcopy(self.p5)
        bad["branches"]["timeout_ungauged"]["full_heading_radius_assigned"] = True
        self.assertTrue(P5.validate(bad))

        bad = deepcopy(self.p5)
        bad["branches"]["timeout_ungauged"]["boundary_cosine_direction_used"] = "LOWER_ENCLOSURE"
        self.assertTrue(P5.validate(bad))

    def test_outer_capture_cannot_be_relabelled_as_final_inner_capture(self):
        self.assertEqual(self.p5["N_outer_words"], 0)
        self.assertIs(self.p5["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"], False)
        bad = deepcopy(self.p5)
        bad["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"] = True
        self.assertTrue(P5.validate(bad))

    def test_438_blockwise_translation_widening_cannot_be_promoted_without_cross_blocks(self):
        translation = {
            "P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": "PASS",
            "modes": {
                "H": {"complete_word_translation_margin_lower": 9.676620313503055e-25},
                "A": {"complete_word_translation_margin_lower": 9.676620313503055e-25},
            },
        }
        post = POST.build(translation)
        self.assertEqual(POST.validate(post), [])
        self.assertIs(post["blockwise_min_is_final_certificate"], False)
        self.assertIs(post["cross_block_budget_is_final_certificate"], False)
        self.assertEqual(post["P4_USABLE_CERTIFICATE_STATUS"], "NOT_ESTABLISHED")
        for mode in ("H", "A"):
            m = post["modes"][mode]
            self.assertTrue(m["cross_block_budget_outward_lower_enclosed"])
            self.assertIs(m["full_state_cross_block_bound_validated"], False)
            self.assertIs(m["full_state_complete_word_cross_blocks_propagated"], False)
            self.assertIs(m["usable_P4_promoted"], False)

        bad = deepcopy(post)
        bad["blockwise_min_is_final_certificate"] = True
        self.assertTrue(POST.validate(bad))

        bad = deepcopy(post)
        bad["cross_block_budget_is_final_certificate"] = True
        self.assertTrue(POST.validate(bad))

    def test_pr441_route_ceiling_must_remain_a_focused_ci_regression(self):
        self.assertIn("tools/ou3_p4_p5_route_ceiling_certificate.py", self.workflow)
        self.assertIn("test_ou3_p4_p5_route_ceiling", self.workflow)

    def test_retired_438_scalar_frontiers_cannot_become_promotion_inputs(self):
        forbidden = (
            "ou3_p4_direct_word_contraction_certificate.py",
            "ou3_p4_nextgen_gain_certificate.py",
            "ou3_p4_nextgen_directional_certificate.py",
            "ou3_p4_nextgen_widened_certificate.py",
            "ou3_p4_thirdgen_combined_certificate.py",
            "ou3_p4_frontier_combined_certificate.py",
        )
        for name in forbidden:
            with self.subTest(name=name):
                self.assertNotIn(name, self.workflow)


if __name__ == "__main__":
    unittest.main()
