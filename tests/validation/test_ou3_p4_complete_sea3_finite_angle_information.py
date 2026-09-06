from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_complete_sea3_finite_angle_information as mod  # noqa: E402


class CompleteSea3FiniteAngleInformationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_complete_sea3_and_actual_rs_are_retained(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertTrue(self.d["component_of_complete_SEA3_full_word"])
        self.assertFalse(self.d["selected_PE_or_four_S_replace_complete_word"])
        self.assertTrue(self.d["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertTrue(
            self.d["actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information"]
        )
        self.assertTrue(self.d["directional_four_S_R_S_regularizer_retained"])

    def test_declared_inner_cells_keep_information_headroom(self):
        self.assertEqual(
            [row["attitude_angle_deg"] for row in self.d["candidate_cells"]],
            [30.0, 25.0, 20.0, 15.0],
        )
        self.assertEqual(self.d["widest_information_cell_deg"], 30.0)
        self.assertTrue(self.d["all_declared_candidate_cells_keep_strict_information"])
        for row in self.d["candidate_cells"]:
            self.assertTrue(row["clears_useful_gate"])
            self.assertGreaterEqual(row["full_H18_information_lambda_min_lower"], 1.0e-18)
        self.assertGreaterEqual(self.d["outer_geometry_angle_rad_retained"], 0.80)
        self.assertTrue(self.d["outer_geometry_sector_separate_from_inner_dissipation_cell"])

    def test_information_headroom_does_not_fake_p4(self):
        self.assertTrue(self.d["magnetometer_radial_state_correction_cancels_exactly"])
        self.assertFalse(self.d["accelerometer_radial_remainder_declared_zero"])
        self.assertTrue(self.d["accelerometer_radial_remainder_requires_signed_Joseph_word_charge"])
        self.assertFalse(self.d["signed_Joseph_complete_word_closed_here"])
        self.assertFalse(self.d["reset_defect_complete_word_closed_here"])
        self.assertFalse(self.d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
