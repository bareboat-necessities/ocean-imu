import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_outer_h_bridge_certificate as BRIDGE


class Ou3P5OuterHBridgeCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = BRIDGE.build()

    def test_bridge_keeps_source_staged_first_S_work(self):
        d = self.d
        self.assertEqual(BRIDGE.validate(d), [])
        self.assertTrue(d["global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate"])
        self.assertTrue(d["finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route"])
        self.assertEqual(d["goLive_covariance_stage"]["status"], "PASS")
        self.assertTrue(d["goLive_covariance_stage"]["S_to_attitude_gain_exact_zero"])
        self.assertEqual(d["first_due_S_gain_stage"]["status"], "PASS")
        p = d["first_due_S_state_prefix_stage"]
        self.assertEqual(p["status"], "PASS_CONDITIONAL")
        self.assertLess(p["first_due_S_induced_attitude_correction_norm_upper_rad"], p["group_helper_limit_rad"])

    def test_bridge_consumes_validated_raw_VR_counterexample(self):
        d = self.d
        a = d["raw_V_R_large_angle_sector_audit"]
        self.assertEqual(a["status"], "DISPROVED_ON_DECLARED_SOURCE_FAMILY")
        self.assertTrue(a["beta_cannot_repair_xi_zero_counterexample"])
        self.assertLess(a["counterexample"]["D_R_deployed_interval"][1], 0.0)
        self.assertTrue(d["raw_V_R_large_angle_sector_retired_as_P5_promotion_route"])
        self.assertTrue(d["source_shaped_Cayley_information_is_primary_full_heading_outer_route"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["raw_V_R_sector_witness_inside_node"])
            self.assertFalse(n["raw_V_R_sector_is_P5_promotion_route"])
            self.assertTrue(n["source_shaped_Cayley_information_outer_sector_required"])
            self.assertEqual(n["source_shaped_Cayley_information_outer_sector_status"], "NOT_ESTABLISHED")

    def test_bridge_uses_full_attitude_gauged_nodes_not_tilt_only_cosines(self):
        d = self.d
        h = d["heading_handoff_contract"]
        self.assertTrue(h["P1_gravity_cosines_are_tilt_only"])
        self.assertGreater(h["gauged_quality_full_cayley_norm_upper"], 0.20)
        self.assertGreater(h["gauged_timeout_full_cayley_norm_upper"], h["gauged_quality_full_cayley_norm_upper"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["P1_gravity_tilt_cosine_not_used_as_full_attitude_cosine"])
            self.assertTrue(n["inside_candidate_outer_cayley_bootstrap"])
            self.assertTrue(n["S_induced_correction_inside_group_helper"])

    def test_yaw_only_quotient_zero_dynamics_is_consumed(self):
        d = self.d
        q = d["yaw_only_quotient_audit"]
        self.assertEqual(q["obstruction_identified"], "PASS")
        self.assertEqual(q["status"], "NOT_ESTABLISHED")
        self.assertTrue(q["witness"]["zero_dynamics_source_word_valid"])
        u = d["ungauged_timeout_route"]
        self.assertTrue(u["yaw_only_quotient_disproved"])
        self.assertFalse(u["full_heading_cayley_bound_available"])
        self.assertIn("AXIAL_GYRO_BIAS", u["required_route"])

    def test_next_obligations_are_the_corrected_metric_and_quotient(self):
        d = self.d
        self.assertEqual(
            d["gauged_full_heading_first_failure"],
            "SOURCE_SHAPED_CAYLEY_INFORMATION_OUTER_SECTOR_NOT_CERTIFIED",
        )
        self.assertEqual(
            d["first_failure"],
            "OBSERVABLE_GRAVITY_ONLY_QUOTIENT_WORD_NOT_CERTIFIED",
        )
        self.assertIn("Cayley/information", d["next_full_heading_numerical_certificate"])
        self.assertIn("axial gyro-bias", d["next_complete_startup_family_certificate"])
        self.assertEqual(d["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIsNone(d["N_H_words"])


if __name__ == "__main__":
    unittest.main()
