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


class OU3CertificateDeadEndNonRegressionTests(unittest.TestCase):
    """Keep already-disproved proof decompositions from returning silently.

    These are not style preferences.  PR #441 established that the scalar
    uniform P4->P5 transport route has a milliradian-scale hard ceiling, while
    PR #438 showed that useful full-word translation widening is only a partial
    block result until the complete full-state cross coupling is enclosed.
    """

    @classmethod
    def setUpClass(cls):
        cls.path = PATH.build()
        cls.p4 = P4.build()
        cls.p5 = P5.build()
        cls.domain = json.loads(DOMAIN.read_text(encoding="utf-8"))

    def test_declared_deployment_domain_cannot_be_shrunk_to_make_certificate_pass(self):
        """Pin the present theorem envelope in the non-shrinking direction."""
        d = self.domain
        self.assertIs(d.get("trajectory_fit"), False)

        s = d["startup"]
        # Startup/P1 admissible disturbance and handoff box: reducing any of
        # these upper bounds would make the theorem easier by excluding states
        # that are currently part of the declared deployment envelope.
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
            "latent_acceleration_error_norm_upper_mps2": 10.0,
            "accelerometer_bias_error_norm_upper_mps2": 0.5,
        }
        for key, baseline in floors.items():
            self.assertGreaterEqual(handoff[key], baseline, key)

        live = d["normal_live"]
        # Raising excitation floors or tightening recurrence is also a domain
        # shrink.  Genuine future widening moves in the opposite direction.
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

    def test_source_path_graph_cannot_be_replaced_by_per_sample_cartesian_corner(self):
        self.assertEqual(PATH.validate(self.path), [])
        self.assertEqual(self.path["P2_SOURCE_PATH_CERTIFICATE"], "PASS")
        self.assertTrue(self.path["path_graph_ready"])
        self.assertGreater(self.path["transition_edges"], 0)
        self.assertGreater(self.path["old_worst_corner_state_count"], 0)
        # The historical weak corner must remain represented, but reachability
        # decides its residence; it may not itself be promoted as a P4 result.
        self.assertFalse(self.path["usable_P4_promoted"])
        self.assertIn("reachable residence/path products", self.path["next_obligation"])
        self.assertNotIn("independently on every sample", self.path["next_obligation"])

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

    def test_ungauged_timeout_cannot_be_given_a_fictitious_full_heading_radius(self):
        u = self.p5["branches"]["timeout_ungauged"]
        self.assertIs(u["full_heading_radius_assigned"], False)
        self.assertEqual(u["attitude_representation"], "GRAVITY_DIRECTION_QUOTIENT")
        bad = deepcopy(self.p5)
        bad["branches"]["timeout_ungauged"]["full_heading_radius_assigned"] = True
        self.assertTrue(P5.validate(bad))

    def test_outer_capture_cannot_be_relabelled_as_final_inner_capture(self):
        self.assertEqual(self.p5["N_outer_words"], 0)
        self.assertIs(self.p5["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"], False)
        bad = deepcopy(self.p5)
        bad["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"] = True
        self.assertTrue(P5.validate(bad))

    def test_438_blockwise_translation_widening_cannot_be_promoted_without_cross_blocks(self):
        # Use a positive synthetic full-word translation result to test the
        # promotion logic without re-running the expensive exact-rational word.
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
            self.assertIs(m["full_state_cross_block_bound_validated"], False)
            self.assertIs(m["full_state_complete_word_cross_blocks_propagated"], False)
            self.assertIs(m["usable_P4_promoted"], False)

        bad = deepcopy(post)
        bad["blockwise_min_is_final_certificate"] = True
        self.assertTrue(POST.validate(bad))

        bad = deepcopy(post)
        bad["cross_block_budget_is_final_certificate"] = True
        self.assertTrue(POST.validate(bad))


if __name__ == "__main__":
    unittest.main()
