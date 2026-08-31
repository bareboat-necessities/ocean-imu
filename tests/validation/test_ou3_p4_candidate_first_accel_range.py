import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_candidate_first_accel_range_v3 as G


class Ou3P4CandidateFirstAccelRangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2, alignment_pieces=16, force_magnitude_pieces=4)

    def test_semantics_validate_and_do_not_change_filter_or_helper_range(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])

    def test_exact_residual_structure_replaces_independent_eta_penalty(self):
        d = self.d
        self.assertTrue(d["candidate_ball_norm_used_without_cartesian_cover_inflation"])
        self.assertTrue(d["rotation_gauge_sets_J_aw_to_identity"])
        self.assertTrue(d["analytic_rank_two_gain_used"])
        self.assertTrue(d["finite_rotation_residual_used_directly"])
        self.assertTrue(d["latent_linear_plus_rotation_cross_combined_before_norm"])
        self.assertFalse(d["independent_accelerometer_eta_penalty_used"])

    def test_H_bound_conservatively_covers_first_A_attitude_gain(self):
        d = self.d
        self.assertTrue(d["H_bias_error_bound_contains_A"])
        self.assertTrue(d["A_first_prefix_attitude_gain_bounded_by_H_gain"])
        self.assertTrue(d["A_structure_proved_before_generic_PSD_boxing"])
        a = d["A_mode_structure"]
        self.assertTrue(a["check_uses_structural_source_matrices_not_psd_tightened_box"])
        self.assertTrue(a["A_bias_seed_cross_exact_zero"])
        self.assertTrue(a["A_bias_transition_cross_exact_zero"])
        self.assertTrue(a["A_bias_process_cross_exact_zero"])
        self.assertTrue(a["A_bias_innovation_addition_isotropic_PSD"])
        self.assertTrue(a["first_prefix_theta_aw_S_to_ba_cross_exact_zero"])
        self.assertTrue(a["A_accelerometer_J_ba_identity"])
        self.assertEqual(a["missing_source_markers"], [])

    def test_tangent_PSD_resolvent_retires_axial_noise_floor(self):
        d = self.d
        self.assertTrue(d["V12D_tangent_PSD_resolvent_used"])
        self.assertFalse(d["PSD_remainder_axial_noise_floor_inverse_used"])

    def test_candidate_ladder_is_complete_and_fail_closed(self):
        d = self.d
        self.assertEqual([r["angle_deg"] for r in d["candidate_rows"]], [30.0, 25.0, 20.0, 15.0])
        for r in d["candidate_rows"]:
            self.assertGreater(r["evaluated_children"], 0)
            self.assertGreater(r["post_prediction_q_upper"], r["candidate_q_upper"])
            self.assertGreaterEqual(r["max_first_accelerometer_correction_norm_upper_rad"], 0.0)
        if d["P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE"] == "PASS":
            self.assertTrue(d["all_candidate_first_accelerometer_ranges_safe"])
            self.assertEqual(d["widest_candidate_first_accel_range_safe_deg"], 30.0)
            self.assertTrue(all(r["max_first_accelerometer_correction_norm_upper_rad"] <= 6.0 for r in d["candidate_rows"]))
        else:
            self.assertTrue(any(r["first_unclosed_child"] is not None for r in d["candidate_rows"]))

    def test_range_stage_does_not_promote_P4(self):
        d = self.d
        self.assertFalse(d["signed_correction_Joseph_reset_propagated_here"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
