from __future__ import annotations

import copy
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))
import ou3_p4_complete_sea3_correction_information_bound as CORR


class CompleteSea3CorrectionInformationBoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CORR.build()

    def test_operation_identity_validates(self):
        self.assertEqual([], CORR.validate(self.d))
        self.assertTrue(CORR.exact_rational_self_test())
        self.assertEqual("d d^T <= J K S K^T <= J P",
                         self.d["full_matrix_correction_inequality"])

    def test_complete_source_and_frozen_p3_are_retained(self):
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", self.d["canonical_source"])
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertTrue(self.d["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertTrue(self.d["same_shipping_P_H_R_K_S_cell_required"])

    def test_no_unproved_numerical_radius_or_floor_is_published(self):
        for key in ("candidate_metric_energy_balls_derived",
                    "source_uniform_covariance_ceiling_certified_here",
                    "posterior_inverse_metric_floor_certified_here",
                    "reset_transport_correction_radius_source_closed",
                    "shipping_Joseph_binding_closed", "P4_promoted_here"):
            self.assertFalse(self.d[key])
            changed = copy.deepcopy(self.d)
            changed[key] = True
            self.assertTrue(CORR.validate(changed), key)
        self.assertNotIn("modes", self.d)
        self.assertTrue(self.d["open_obligations"])

    def test_attitude_marginal_is_not_inverse_metric_bound(self):
        # P=[[1,r],[r,1]], r=9/10: P_00=1, but (P^-1)_00=100/19>1.
        correlation = F(9, 10)
        inverse_attitude = 1 / (1 - correlation * correlation)
        self.assertGreater(inverse_attitude, 1)
        self.assertFalse(self.d["posterior_inverse_metric_floor_certified_here"])

    def test_reset_isometry_does_not_bound_euclidean_covariance(self):
        # G=I+[d]_x/2 is an exact metric congruence, not an orthogonal map.
        d = F(1, 10)
        reset_transverse_variance = 1 + d*d/4
        self.assertGreater(reset_transverse_variance, 1)
        self.assertFalse(self.d["source_uniform_covariance_ceiling_certified_here"])

    def test_candidate_geometry_is_retained_without_packet_or_tube_shortcut(self):
        self.assertEqual([30.0, 25.0, 20.0, 15.0], self.d["candidate_angles_deg"])
        self.assertFalse(self.d["old_scalar_Riccati_tube_margin_consumed"])
        self.assertFalse(self.d["packet_count_multiplier_used"])
        self.assertFalse(self.d["independent_global_correction_radius_assumed"])
        self.assertFalse(self.d["complete_word_nonlinear_dissipation_closed_here"])


if __name__ == "__main__":
    unittest.main()
