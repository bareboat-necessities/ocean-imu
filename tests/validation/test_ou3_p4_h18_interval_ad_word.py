from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_interval_ad as AD
import ou3_p4_h18_interval_ad_word as H18


class H18IntervalADWordTests(unittest.TestCase):
    def test_exact_accelerometer_residual_has_expected_local_H_blocks(self):
        z = [AD.independent(Interval.point(0.0), i, H18.N) for i in range(H18.N)]
        f = [Interval.point(0.0), Interval.point(0.0), Interval.point(7.0)]
        r = H18._rotation_residual_acc(z, f)
        J = AD.jacobian(r)
        H = H18._H_acc_canonical(f)
        for i in range(3):
            for j in list(range(0, 3)) + list(range(15, 18)):
                self.assertLessEqual(J[i][j].lo, H[i][j].lo)
                self.assertGreaterEqual(J[i][j].hi, H[i][j].hi)

    def test_exact_magnetometer_residual_has_expected_local_attitude_H(self):
        z = [AD.independent(Interval.point(0.0), i, H18.N) for i in range(H18.N)]
        m = [Interval.point(2.0), Interval.point(0.0), Interval.point(9.0)]
        r = H18._rotation_residual_mag(z, m)
        J = AD.jacobian(r)
        H = H18._H_mag_canonical(m)
        for i in range(3):
            for j in range(3):
                self.assertLessEqual(J[i][j].lo, H[i][j].lo)
                self.assertGreaterEqual(J[i][j].hi, H[i][j].hi)

    def test_S_residual_is_exact_negative_identity_on_S_block(self):
        z = [AD.independent(Interval.point(0.0), i, H18.N) for i in range(H18.N)]
        J = AD.jacobian(H18._residual_S(z))
        for i in range(3):
            for j in range(H18.N):
                expected = -1.0 if j == 12 + i else 0.0
                self.assertEqual(J[i][j].lo, expected)
                self.assertEqual(J[i][j].hi, expected)

    def test_conditioning_similarity_preserves_identity(self):
        Z = Interval.point(0.0)
        O = Interval.point(1.0)
        J = [[O if i == j else Z for j in range(H18.N)] for i in range(H18.N)]
        scale2 = [float(i + 1) for i in range(H18.N)]
        C = H18._conditioned_jacobian(J, scale2)
        for i in range(H18.N):
            for j in range(H18.N):
                expected = 1.0 if i == j else 0.0
                self.assertLessEqual(C[i][j].lo, expected)
                self.assertGreaterEqual(C[i][j].hi, expected)

    def _gamma_fixture(self):
        return {
            "schema": H18.SCHEMA,
            "qualification": "OU3_P4_H18_INTERVAL_AD_WORD_SCREEN",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "filter_changed": False,
            "dimension": 18,
            "outer_angle_rad": 0.80,
            "outer_ball_box_cover_total": 1,
            "interval_AD_used_for_state_return_map": True,
            "finite_difference_used": False,
            "deployed_quaternion_generalized_jacobian_used": True,
            "full_18_state_cross_derivatives_retained": True,
            "canonical_common_rotation_PE_gauge_used": True,
            "all_source_vector_orientation_covariance_correlations_checked": False,
            "mandatory_four_S_events_used": True,
            "mandatory_two_packet_vector_PE_used": True,
            "optional_accepted_branch_family_between_required_events_checked": False,
            "aw_covariance_sync_overapproximated_at_every_prefix": True,
            "actual_per_node_Sigma_KF_whitening_used": False,
            "P3_computational_congruence_used_for_screening_only": True,
            "P3_delta_used_as_nonlinear_radius": False,
            "source_graph_all_reachable_edges_checked": False,
            "H18_COMPLETE_SOURCE_EDGE_CONTRACTION_ESTABLISHED_HERE": False,
            "P4_USABLE_CERTIFICATE_PROMOTED": False,
            "H18_SCREEN_GAMMA_LT_ONE": True,
            "full_word_horizon_checked": True,
            "all_outer_ball_cells_checked": True,
            "max_endpoint_P3_congruence_conditioned_norm_upper": 0.9,
            "max_prefix_P3_congruence_conditioned_norm_upper": 1.2,
            "first_failure": None,
            "failures": [],
        }

    def test_screen_validator_is_fail_closed_against_premature_promotion(self):
        d = self._gamma_fixture()
        self.assertEqual(H18.validate(d), [])
        d["P4_USABLE_CERTIFICATE_PROMOTED"] = True
        self.assertIn("P4_USABLE_CERTIFICATE_PROMOTED is not false", H18.validate(d))

    def test_subunity_gamma_requires_complete_word_ball_and_no_failure(self):
        d = self._gamma_fixture()
        d["full_word_horizon_checked"] = False
        self.assertIn(
            "sub-unity H18 gamma claimed without complete word horizon",
            H18.validate(d),
        )

        d = self._gamma_fixture()
        d["all_outer_ball_cells_checked"] = False
        self.assertIn(
            "sub-unity H18 gamma claimed without complete outer-ball coverage",
            H18.validate(d),
        )

        d = self._gamma_fixture()
        d["first_failure"] = {"error": "synthetic numerical failure"}
        self.assertIn(
            "sub-unity H18 gamma claimed despite numerical failure",
            H18.validate(d),
        )

    def test_truncated_schedule_reports_actual_mandatory_event_counts(self):
        # The schedule itself is cheap to inspect and demonstrates why a
        # truncated smoke must not claim the complete four-S/two-vector family.
        words = H18.WORDS.build(H18.DEFAULT_DOMAIN)
        sched = H18._mandatory_schedule(words, samples=1, h=0.005)
        self.assertEqual(len(sched["S_steps"]), 1)
        self.assertEqual(len(sched["vector_steps"]), 1)


if __name__ == "__main__":
    unittest.main()
