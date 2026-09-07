#!/usr/bin/env python3
from fractions import Fraction as F
import sys
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
STABILITY = REPO / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p4_complete_word_endpoint_transport as ENDPOINT
import ou3_sea3_riccati_metric_p4 as P4


def mat(rows):
    return [[F(x) for x in row] for row in rows]


def vec(xs):
    return [F(x) for x in xs]


class CompleteWordEndpointTransportTests(unittest.TestCase):
    def test_full_shift_telescope_leaves_only_accelerometer_interior_block(self):
        I2 = mat([[1, 0], [0, 1]])
        events = [
            {
                "kind": "prediction",
                "C": mat([[1, 1], [0, 1]]),
                "L": mat([[1, 1], [0, 1]]),
                "rho": vec([F(1, 10), 0]),
            },
            {
                "kind": "S_zero",
                "C": mat([[F(1, 2), 0], [0, 1]]),
                "L": I2,
                "rho": vec([0, F(1, 20)]),
            },
            {
                "kind": "accelerometer",
                "C": mat([[1, 0], [0, F(1, 2)]]),
                "L": I2,
                "rho": vec([F(-1, 30), 0]),
            },
            {
                "kind": "hybrid",
                "C": mat([[1, 0], [0, 1], [0, 0]]),
                "L": mat([[1, 0], [0, 1], [0, 0]]),
                "rho": vec([0, 0, F(1, 40)]),
            },
            {
                "kind": "prediction",
                "C": mat([[1, 0, 1], [0, 1, 0], [0, 0, 1]]),
                "L": mat([[1, 0, 1], [0, 1, 0], [0, 0, 1]]),
                "rho": vec([0, F(-1, 50), 0]),
            },
        ]
        embeddings = [
            mat([[0], [1]]),
            mat([[0], [1]]),
            mat([[0], [1]]),
            mat([[0], [1]]),
            mat([[0], [1], [0]]),
            mat([[0], [1], [0]]),
        ]
        eps = [[F(1, 7)], [F(2, 7)], [F(-1, 9)], [F(1, 11)], [F(3, 13)], [F(-2, 15)]]

        d = ENDPOINT.endpoint_decomposition(events, embeddings, eps)
        self.assertEqual(d["identity_residual"], vec([0, 0, 0]))
        self.assertEqual(d["direct_defect"], d["decomposed_defect"])
        self.assertEqual(d["interior_event_indices"], [2])
        self.assertEqual(len(d["accelerometer_suffix_blocks"]), 1)

    def test_rectangular_hybrid_cancels_when_C_equals_L(self):
        events = [{
            "kind": "hybrid",
            "C": mat([[1, 0], [0, 1], [1, -1]]),
            "L": mat([[1, 0], [0, 1], [1, -1]]),
            "rho": vec([F(1, 5), 0, F(-1, 7)]),
        }]
        embeddings = [mat([[0], [1]]), mat([[0], [1], [-1]])]
        eps = [[F(2, 9)], [F(-1, 4)]]
        d = ENDPOINT.endpoint_decomposition(events, embeddings, eps)
        self.assertEqual(d["identity_residual"], vec([0, 0, 0]))
        self.assertEqual(d["interior_event_indices"], [])

    def test_nonaccelerometer_nonzero_interior_transport_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "non-accelerometer"):
            ENDPOINT.endpoint_decomposition(
                [{
                    "kind": "S_zero",
                    "C": mat([[1, 0], [0, F(1, 2)]]),
                    "L": mat([[1, 0], [0, 1]]),
                    "rho": vec([0, 0]),
                }],
                [mat([[0], [1]]), mat([[0], [1]])],
                [[F(1, 3)], [F(1, 4)]],
            )

    def test_master_energy_identity_is_exact_but_nonpromoting(self):
        J0 = mat([[2, F(1, 3)], [F(1, 3), 3]])
        JN = mat([[4, F(1, 5)], [F(1, 5), 5]])
        M = mat([[F(1, 2), F(1, 7)], [0, F(1, 3)]])
        phi0 = vec([F(2, 5), F(-3, 7)])
        defect = vec([F(1, 11), F(-1, 13)])
        d = ENDPOINT.master_energy_identity(J0, JN, M, phi0, defect)
        self.assertEqual(d["identity_residual"], F(0))
        self.assertEqual(d["direct_delta"], d["total_delta"])
        self.assertEqual(
            d["total_delta"],
            -d["linear_decrease"] + d["cross_term"] + d["defect_energy"],
        )

    def test_endpoint_status_retains_complete_sea3_and_forbids_shortcuts(self):
        d = ENDPOINT.build()
        self.assertEqual(ENDPOINT.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_consumed"], 1.0e-18)
        self.assertTrue(d["actual_applied_per_axis_RS_retained_in_word"])
        self.assertTrue(d["accepted_accelerometer_is_only_interior_epsilon_event_class"])
        self.assertFalse(d["packetwise_norm_sum_used"])
        self.assertFalse(d["packet_count_multiplier_used"])
        self.assertFalse(d["correction_radius_claim_used"])
        self.assertFalse(d["inverse_metric_floor_claim_used"])
        self.assertFalse(d["P4_promoted_here"])

    def test_endpoint_transport_is_subordinate_to_canonical_finite_state_p4(self):
        d = P4.build()
        self.assertEqual(P4.validate(d), [])
        self.assertEqual(
            d["canonical_P4_architecture"],
            "FINITE_STATE_COMPLETE_SEA3_QUADRATIC_ENDPOINT_AND_PREFIX",
        )
        self.assertEqual(d["paper_Lyapunov_function"], "V(e,zeta)=e^T M(zeta)e")
        self.assertTrue(d["finite_state_endpoint_dissipation_required"])
        self.assertTrue(d["finite_state_prefix_gain_required"])
        self.assertFalse(d["differential_pullback_used_as_replacement_P4"])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        self.assertFalse(d["P5_MAY_START"])


if __name__ == "__main__":
    unittest.main()
