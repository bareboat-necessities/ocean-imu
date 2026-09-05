from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_full_word_riccati_backend as mod  # noqa: E402
from ou3_interval import matrix_point  # noqa: E402


class Sea3FullWordRiccatiBackendTest(unittest.TestCase):
    def test_shipping_source_parity(self):
        p = mod.shipping_source_parity()
        self.assertTrue(p)
        self.assertTrue(all(p.values()), p)

    def test_prediction_measurement_floor_preserve_exact_decomposition(self):
        state = mod.initialize(matrix_point([
            [2.0, 0.1],
            [0.1, 3.0],
        ]))
        self.assertTrue(mod.decomposition_identity_enclosed(state))

        mod.predict(
            state,
            matrix_point([[1.0, 0.1], [0.0, 1.0]]),
            matrix_point([[0.02, 0.0], [0.0, 0.01]]),
        )
        self.assertTrue(mod.decomposition_identity_enclosed(state))

        mod.joseph_measurement(
            state,
            matrix_point([[1.0, 0.0]]),
            matrix_point([[0.5]]),
        )
        self.assertTrue(mod.decomposition_identity_enclosed(state))

        psi_before = [[x for x in row] for row in state.Psi]
        mod.add_psd_floor(
            state,
            matrix_point([[0.0, 0.0], [0.0, 0.03]]),
        )
        self.assertTrue(mod.decomposition_identity_enclosed(state))
        self.assertEqual(state.Psi, psi_before)

    def test_exact_M_delta_identities_cover_every_shipping_event_type(self):
        ids = mod.contraction_preservation_identities()
        self.assertEqual(ids["definition"], "M_delta = Omega - delta P")
        self.assertEqual(
            ids["prediction"],
            "M_delta' = F M_delta F^T + (1-delta) Q",
        )
        self.assertEqual(
            ids["Joseph"],
            "M_delta' = A M_delta A^T + (1-delta) K R K^T",
        )
        self.assertEqual(
            ids["PSD_floor"],
            "M_delta' = M_delta + (1-delta) Delta",
        )
        self.assertTrue(ids["positive_semidefinite_inputs_preserve_M_delta"])
        self.assertTrue(ids["is_exact_algebra_inside_joint_P_Psi_Omega_object"])
        self.assertFalse(ids["may_omit_events_needed_to_first_establish_delta"])

    def test_M_delta_event_images_equal_direct_joint_recursion(self):
        t = mod._preservation_self_test()
        self.assertTrue(t["prediction_identity_enclosed"])
        self.assertTrue(t["Joseph_identity_enclosed"])
        self.assertTrue(t["PSD_floor_identity_enclosed"])
        self.assertTrue(t["all_event_identities_enclosed"])
        self.assertTrue(t["kernel_identity_test_only_not_P3"])

    def test_only_full_matrix_omega_minus_delta_p_gate_is_exposed(self):
        state = mod.initialize(matrix_point([
            [2.0, 0.0],
            [0.0, 3.0],
        ]))
        mod.predict(
            state,
            matrix_point([[1.0, 0.0], [0.0, 1.0]]),
            matrix_point([[0.1, 0.0], [0.0, 0.1]]),
        )
        cert = mod.certify_contraction(state, mod.USEFUL_GATE)
        self.assertEqual(cert["dimension"], 2)
        self.assertEqual(cert["delta"], 1e-18)
        self.assertTrue(cert["decomposition_identity_enclosed"])
        self.assertTrue(cert["Omega_minus_delta_P_ldlt_closed"])
        self.assertTrue(cert["pass"])

    def test_backend_validation_does_not_promote_ou3_p3(self):
        self.assertEqual(mod.validate_backend(), [])
        test = mod._self_test()
        self.assertTrue(test["kernel_self_test_only_not_P3"])
        self.assertTrue(test["shipping_source_parity_pass"])
        self.assertTrue(test["contraction_preservation"]["all_event_identities_enclosed"])


if __name__ == "__main__":
    unittest.main()
