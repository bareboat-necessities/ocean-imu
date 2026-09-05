from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_complete_source as source  # noqa: E402
import ou3_sea3_full_normal_live_word as word  # noqa: E402
from ou3_interval import Interval  # noqa: E402


class Sea3RSFullWordContractTest(unittest.TestCase):
    def test_actual_spectral_mse_rs_is_mandatory(self):
        d = source.build()
        self.assertEqual(source.validate(d), [])
        rs = d["R_S_regularizer"]
        self.assertEqual(rs["deployed_law"], "SpectralMSE")
        self.assertEqual(rs["axis_std_factors"], [0.72, 0.72, 1.0])
        self.assertTrue(rs["actual_applied_R_S_required_at_every_due_S_update"])
        self.assertTrue(rs["all_due_S_updates_remain_in_full_word"])
        self.assertTrue(rs["full_P_column_S_cross_covariance_action_required"])

    def test_literal_rs_matrix_squares_supplied_applied_standard_deviation(self):
        r = [Interval.point(7.2), Interval.point(7.2), Interval.point(10.0)]
        R = word.R_S_zero(r)
        self.assertEqual(R[0][0], Interval.point(7.2).square())
        self.assertEqual(R[1][1], Interval.point(7.2).square())
        self.assertEqual(R[2][2], Interval.point(10.0).square())


if __name__ == "__main__":
    unittest.main()
