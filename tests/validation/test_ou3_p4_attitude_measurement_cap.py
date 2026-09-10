#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_attitude_measurement_cap as CAP


class AttitudeMeasurementCapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CAP.build()

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(CAP.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])

    def test_prior_independent_cap_is_recorded_as_refuted(self):
        self.assertFalse(self.d["prior_independent_transverse_cap_holds_for_the_deployed_filter"])
        ref = self.d["refutation"]
        self.assertFalse(ref["claim_holds"])
        self.assertGreater(ref["achieved_transverse_attitude_marginal"],
                           ref["claimed_prior_independent_cap"])
        self.assertGreater(ref["exceedance_factor"], 1e6)

    def test_refutation_needs_the_latent_block_to_bite(self):
        # With attitude alone in the residual the scalar identity is correct, so
        # the refutation must come from the shared residual and not from a bug.
        alone = CAP.refutation(latent_prior=0.0)
        self.assertTrue(alone["claim_holds"])
        self.assertFalse(alone["refutation_is_outside_the_noise_floor"])
        self.assertLessEqual(alone["achieved_transverse_attitude_marginal"],
                             alone["claimed_prior_independent_cap"]
                             + alone["cancellation_noise_floor"])
        # And it must grow monotonically with the latent prior.
        prev = 0.0
        for lp in (1.0, 1.0e2, 1.0e4, 1.0e8):
            got = CAP.refutation(latent_prior=lp)["achieved_transverse_attitude_marginal"]
            self.assertGreater(got, prev, lp)
            prev = got

    def test_conditional_cap_dominates_the_deployed_posterior(self):
        # Direct check of the variational bound on the deployed nine-state
        # residual structure H = [-skew(f), R^T, I]: build a posterior by the
        # Joseph/innovation form and confirm the transverse marginal obeys it.
        import math
        import random

        rng = random.Random(20260910)
        f_norm, sigma = 9.80665, 0.2
        f = [0.0, 0.0, -f_norm]
        Hth = [[-x for x in row] for row in CAP._skew(f)]
        worst = 0.0
        for _ in range(400):
            # (theta, a_w) with the latent block at identity rotation.
            H = [Hth[i] + [1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
            P = [[0.0] * 6 for _ in range(6)]
            lam = 0.0
            for i in range(6):
                v = 10.0 ** rng.uniform(-3, 4)
                P[i][i] = v
                if i >= 3:
                    lam = max(lam, v)
            HP = CAP._matmul(H, P)
            S = CAP._matmul(HP, CAP._transpose(H))
            for i in range(3):
                S[i][i] += sigma * sigma
            cols = [CAP._solve(S, [HP[i][j] for i in range(3)]) for j in range(6)]
            Pp = [[P[i][j] - sum(HP[k][i] * cols[j][k] for k in range(3)) for j in range(6)]
                  for i in range(6)]
            bound = CAP.conditional_cap(sigma, f_norm, lam)
            for axis in (0, 1):  # both directions transverse to f
                self.assertLessEqual(Pp[axis][axis], bound * (1.0 + 1e-9))
                worst = max(worst, Pp[axis][axis] / bound)
        self.assertGreater(worst, 0.1, "sweep never approached the bound")
        self.assertTrue(math.isfinite(worst))

    def test_conditional_cap_is_monotone_and_matches_its_closed_form(self):
        a = CAP.conditional_cap(0.2, 9.8, 1.0)
        b = CAP.conditional_cap(0.2, 9.8, 2.0)
        self.assertGreater(b, a)
        self.assertGreater(CAP.conditional_cap(0.4, 9.8, 1.0), a)
        self.assertGreater(a, CAP.conditional_cap(0.2, 12.0, 1.0))
        self.assertAlmostEqual(a, (0.2 * 0.2 + 2.0 * 1.0) / (9.8 * 9.8), places=12)
        for bad in ((0.0, 9.8, 1.0), (0.2, 0.0, 1.0), (0.2, 9.8, -1.0)):
            with self.assertRaises(ValueError):
                CAP.conditional_cap(*bad)

    def test_sweep_finds_no_violation_and_shows_the_bound_is_tight(self):
        sw = self.d["sweep"]
        self.assertEqual(sw["violations"], 0)
        self.assertGreater(sw["numerically_valid_cases_kept"], sw["cases"] // 2)
        self.assertGreater(sw["ill_conditioned_cases_rejected"], 0,
                           "no rejections means the validity gate is not doing anything")
        self.assertLessEqual(sw["max_achieved_over_bound"], 1.0)
        self.assertGreaterEqual(sw["max_achieved_over_bound"],
                                CAP.CONDITIONAL_BOUND_TIGHTNESS_FLOOR)
        self.assertIn("R^T", sw["residual_structure"])
        # The worst case must be a real one, with its own numbers, so the
        # tightness figure can be rechecked by hand.
        wc = sw["worst_case"]
        self.assertAlmostEqual(
            wc["bound"],
            CAP.conditional_cap(wc["sigma_accelerometer"], wc["specific_force_norm"],
                                wc["joint_latent_bias_lambda_max"]),
            places=15)
        self.assertLessEqual(wc["achieved_transverse_marginal"], wc["bound"])

    def test_sweep_is_deterministic_and_gated_on_the_variational_identity(self):
        a = CAP.sweep(400, seed=7)
        b = CAP.sweep(400, seed=7)
        self.assertEqual(a, b)
        self.assertNotEqual(CAP.sweep(400, seed=8)["worst_case"], a["worst_case"])
        self.assertEqual(a["numerically_valid_cases_kept"] + a["ill_conditioned_cases_rejected"],
                         a["cases"])
        self.assertEqual(a["violations"], 0)

    def test_cap_does_not_close_the_correction_domain(self):
        self.assertFalse(self.d["cap_closes_the_correction_reset_domain"])
        self.assertFalse(self.d["envelope_shown_to_be_loose_by_this_argument"])
        self.assertTrue(self.d["one_shot_information_route_is_a_dead_end_for_the_correction_domain"])
        self.assertTrue(self.d["attitude_covariance_bound_needs_uniform_observability_over_a_window"])
        self.assertIn("uniform observability", self.d["next_obligation"])

    def test_lambda_max_comes_from_the_certified_ceiling_not_an_assumption(self):
        self.assertEqual(self.d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertGreater(self.d["joint_latent_bias_lambda_max"], 0.0)
        self.assertGreaterEqual(self.d["joint_latent_bias_lambda_max"],
                                self.d["latent_lambda_max"])
        self.assertGreaterEqual(self.d["joint_latent_bias_lambda_max"],
                                self.d["bias_lambda_max"])

    def test_validate_rejects_a_resurrected_prior_independent_claim(self):
        d = dict(self.d)
        d["prior_independent_transverse_cap_holds_for_the_deployed_filter"] = True
        self.assertIn(
            "prior_independent_transverse_cap_holds_for_the_deployed_filter not false",
            CAP.validate(d))
        d = dict(self.d)
        d["refutation"] = dict(self.d["refutation"], claim_holds=True)
        self.assertIn("refutation does not refute the claimed cap", CAP.validate(d))
        d = dict(self.d)
        d["conditional_transverse_cap_rad2"] = 1.0e-4
        self.assertIn("conditional cap is suspiciously small; recheck before consuming it",
                      CAP.validate(d))


if __name__ == "__main__":
    unittest.main()
