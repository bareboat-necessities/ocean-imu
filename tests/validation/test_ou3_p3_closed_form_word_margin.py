"""Contracts for the closed-form whole-word translation margin."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_closed_form_word_margin as M  # noqa: E402
import ou3_source_reachable_matrix_p3 as BASE  # noqa: E402


class ShippingAgreement(unittest.TestCase):
    def test_cadence_matches_the_repo_source_parser(self):
        k = M._consts()
        sched = BASE.source_schedule()
        for tau in (0.4, 1.1, 3.0, 6.0, 12.0):
            want = min(max(sched["pseudo_ratio"] * tau,
                           sched["pseudo_min_s"]), sched["pseudo_max_s"])
            self.assertEqual(M.cadence_s(tau, k), want)

    def test_R_S_is_the_clamped_deployed_law_not_a_free_cell_range(self):
        k = M._consts()
        # The law derives R_S from tau and sigma; it must respect both clamps.
        for tau, sigma in ((0.4, 0.05), (3.0, 1.0), (12.0, 4.0), (12.0, 6.0)):
            rs = M.rs_target(tau, sigma, k)
            self.assertGreaterEqual(rs, k["min_rs"])
            self.assertLessEqual(rs, k["max_rs"])
        # Monotone increasing in sigma and in tau, as the closed form requires.
        self.assertLess(M.rs_target(3.0, 0.05, k), M.rs_target(3.0, 1.0, k))
        self.assertLess(M.rs_target(1.1, 1.0, k), M.rs_target(6.0, 1.0, k))

    def test_initial_covariance_is_the_shipping_constructor_value(self):
        self.assertEqual((M.SIGMA_V0, M.SIGMA_P0, M.SIGMA_S0), (1.0, 20.0, 50.0))


class ClosedFormStructure(unittest.TestCase):
    def test_gramian_grows_with_more_observations(self):
        k = M._consts()
        # A faster cadence means more S observations, hence strictly more
        # information in the S direction.
        fast, n_fast = M._gramian(1.0, 1.0, 10.0, k, inflate=False)
        slow, n_slow = M._gramian(12.0, 1.0, 10.0, k, inflate=False)
        self.assertGreater(n_fast, n_slow)
        self.assertGreater(fast[2][2].lo, slow[2][2].lo)

    def test_uninflated_gramian_carries_more_information_than_inflated(self):
        k = M._consts()
        bare, _ = M._gramian(6.0, 1.0, 10.0, k, inflate=False)
        infl, _ = M._gramian(6.0, 1.0, 10.0, k, inflate=True)
        # The floor credits full-strength observations, so its Gramian dominates.
        for i in range(3):
            self.assertGreaterEqual(bare[i][i].lo, infl[i][i].lo)

    def test_floor_never_exceeds_ceiling_on_any_channel(self):
        k = M._consts()
        import ou3_p4_source_node_cells as NODES
        for node in NODES.build()["nodes"][::173]:
            r = M.node_margin(node, k)
            self.assertTrue(r["validated"], r.get("reason"))
            for x in r["channel_ratios"]:
                self.assertGreater(x, 0.0)
                self.assertLessEqual(x, 1.0 + 1e-9)

    def test_inversion_is_validated_not_assumed(self):
        # A singular Gramian must be reported, never silently inverted.
        G = [[M._I(0.0) for _ in range(3)] for _ in range(3)]
        self.assertIsNone(M._inv3(G))


class NonPromotion(unittest.TestCase):
    def test_artifact_declares_what_it_is_and_is_not(self):
        d = M.build(stride=401)
        self.assertTrue(d["non_promoting"])
        self.assertFalse(d["certifies_theorem_stage"])
        self.assertTrue(d["interval_certified"])
        self.assertFalse(d["stepwise_riccati_recursion_used"])
        self.assertEqual(d["canonical_gate"], 1e-18)
        self.assertEqual(M.validate(d), [])

    def test_validate_rejects_a_tampered_gate_or_promotion_flag(self):
        d = M.build(stride=401)
        for key, bad in (("canonical_gate", 1e-12),
                         ("certifies_theorem_stage", True),
                         ("non_promoting", False),
                         ("stepwise_riccati_recursion_used", True)):
            t = dict(d)
            t[key] = bad
            self.assertTrue(M.validate(t), f"tampering with {key} was not caught")


if __name__ == "__main__":
    unittest.main()
