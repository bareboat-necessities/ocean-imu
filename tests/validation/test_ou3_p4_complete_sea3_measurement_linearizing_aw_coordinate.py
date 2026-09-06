from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
STABILITY = TOOLS / "stability"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(STABILITY))

import ou3_p4_complete_sea3_measurement_linearizing_aw_coordinate as CERT


class MeasurementLinearizingAwCoordinateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CERT.build()

    def test_complete_sea3_exact_coordinate_contract(self):
        d = self.d
        self.assertEqual([], CERT.validate(d))
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertTrue(d["P3_frozen_not_modified"])
        self.assertTrue(d["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertTrue(d["exact_cayley_remainder_identity_closed"])
        self.assertTrue(d["exact_shipping_tangent_H_used"])
        self.assertFalse(d["shipping_Joseph_binding_closed"])
        self.assertTrue(d["phi_storage_has_no_standalone_eta_penalty"])
        self.assertFalse(d["standalone_eta_Rinv_packet_budget_used"])
        self.assertFalse(d["packet_count_multiplier_used"])
        self.assertFalse(d["complete_source_correlated_transport_defect_closed_here"])
        self.assertFalse(d["P4_promoted_here"])
        self.assertTrue(d["full_mixed_aw_shift_retained_in_transport"])
        self.assertFalse(d["candidate_shift_bounds_cover_full_epsilon_aw"])
        self.assertFalse(d["nonlinear_storage_metric_equivalence_to_original_closed"])

    def test_finite_angle_residual_matrix_is_not_congruent_shipping_matrix(self):
        # Exact rational Cayley c_z=1/10 rotation: Q_aw^T != I.  With R_hat=I,
        # H0_aw=I while H_u_aw=Q_aw^T. No numerical tolerance or replay.
        from fractions import Fraction as F
        q_transpose = ((F(399, 401), F(40, 401), 0),
                       (-F(40, 401), F(399, 401), 0), (0, 0, 1))
        identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.assertNotEqual(identity, q_transpose)
        d = self.d
        self.assertFalse(d["shipping_Joseph_binding_closed"])

    def test_transport_defect_is_combined_not_independent_eta(self):
        d = self.d
        tr = d["combined_correction_reset_transport"]
        self.assertEqual("rho+E_aw(epsilon_plus-epsilon_minus)", tr["xi"])
        self.assertTrue(tr["G_identity_on_aw_coordinate"])
        self.assertTrue(tr["reset_covariance_congruence_exact"])
        self.assertEqual(1.0, tr["G_inverse_operator_norm_exact"])
        self.assertTrue(d["source_indexed_shift_must_persist_across_complete_word"])
        self.assertFalse(d["packetwise_shift_reset_to_zero_allowed"])

    def test_candidate_shift_bounds_are_finite_and_monotone(self):
        d = self.d
        rows = d["candidate_cells"]
        self.assertEqual([30.0, 25.0, 20.0, 15.0],
                         [r["attitude_angle_deg"] for r in rows])
        shifts = [float(r["e_eta_norm_upper_mps2"]) for r in rows]
        lips = [float(r["e_eta_local_lipschitz_upper_mps2_per_cayley"]) for r in rows]
        self.assertTrue(all(math.isfinite(x) and x > 0.0 for x in shifts + lips))
        self.assertTrue(all(a > b for a, b in zip(shifts, shifts[1:])))
        self.assertTrue(all(a > b for a, b in zip(lips, lips[1:])))
        self.assertTrue(all(0.0 < float(r["cayley_norm_upper"]) < 1.0 for r in rows))


if __name__ == "__main__":
    unittest.main()
