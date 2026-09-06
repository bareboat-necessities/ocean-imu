from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
sys.path.insert(0, str(STABILITY))

import ou3_p4_complete_sea3_correction_information_bound as CORR


class CompleteSea3CorrectionInformationBoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CORR.build()

    def test_certificate_validates(self):
        self.assertEqual([], CORR.validate(self.d))

    def test_complete_source_and_frozen_p3_are_retained(self):
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", self.d["canonical_source"])
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertTrue(self.d["P3_conditional_complete_SEA3_consumed"])
        self.assertTrue(self.d["complete_SEA3_word_retained"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertTrue(self.d["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertTrue(self.d["actual_RS_enters_same_Riccati_tube"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["source_family_replaced"])

    def test_correction_radius_is_derived_from_same_operation_information(self):
        self.assertTrue(self.d["same_shipping_P_H_R_K_S_cell_used"])
        self.assertTrue(self.d["same_operation_correction_radius_derived_from_information"])
        self.assertEqual(
            "H^T S^-1 H <= P^-1",
            self.d["measurement_information_operator_inequality"],
        )
        self.assertIn("E_theta K y", self.d["attitude_correction_information_inequality"])
        self.assertFalse(self.d["independent_global_correction_radius_assumed"])

    def test_each_candidate_has_positive_derived_energy_ball_and_posterior_floor(self):
        self.assertEqual([30.0, 25.0, 20.0, 15.0], self.d["candidate_angles_deg"])
        for mode in ("H", "A"):
            rows = self.d["modes"][mode]
            self.assertEqual(4, len(rows))
            previous_q = math.inf
            for row in rows:
                ptheta = float(row["attitude_covariance_lambda_max_upper"])
                floor = float(row["post_measurement_attitude_covariance_floor"])
                nu = float(row["derived_metric_energy_radius_upper"])
                q = float(row["candidate_cayley_norm_upper"])
                self.assertTrue(math.isfinite(ptheta) and ptheta > 0.0)
                self.assertTrue(math.isfinite(floor) and floor > 0.0)
                self.assertTrue(math.isfinite(nu) and nu > 0.0)
                self.assertLess(q, previous_q)
                previous_q = q
                self.assertAlmostEqual(q * q / ptheta, nu, delta=1e-12 * max(1.0, nu))

    def test_reset_cost_uses_posterior_floor_not_independent_packet_budget(self):
        self.assertTrue(self.d["complete_SEA3_Riccati_tube_supplies_post_measurement_floor"])
        for mode in ("H", "A"):
            for row in self.d["modes"][mode]:
                diag = row["reset_endpoint_diagnostic"]
                self.assertGreater(diag["reset_attitude_defect_norm_upper"], 0.0)
                self.assertGreater(diag["reset_defect_metric_cost_upper"], 0.0)
                self.assertFalse(diag["additive_packet_budget_used"])
        self.assertFalse(self.d["packet_count_multiplier_used"])
        self.assertFalse(self.d["standalone_eta_Rinv_budget_used"])
        self.assertFalse(self.d["endpoint_source_word_enumeration_used"])

    def test_helper_closes_only_reset_radius_source_not_p4(self):
        self.assertTrue(self.d["reset_transport_correction_radius_source_closed"])
        self.assertFalse(self.d["source_indexed_e_eta_transition_closed_here"])
        self.assertFalse(self.d["complete_word_nonlinear_dissipation_closed_here"])
        self.assertFalse(self.d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
