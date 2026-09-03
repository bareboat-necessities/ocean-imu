from pathlib import Path
import math
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p3_whole_word_lower_feasibility as F


class SchedulerSemantics(unittest.TestCase):
    """The diagnostic is only meaningful if it reproduces the shipping scheduler."""

    def test_cadence_is_tau_scaled_and_clamped(self):
        self.assertAlmostEqual(F.cadence_s(1.1), 0.015, places=6)
        # Both safety clamps are reachable from the declared tau envelope ends.
        self.assertEqual(F.cadence_s(1.0e-6), F.PSEUDO_MIN_S)
        self.assertEqual(F.cadence_s(1.0e6), F.PSEUDO_MAX_S)

    def test_periodic_update_due_retains_remainder_and_never_fires_at_zero(self):
        period, elapsed, fires = 0.02, 0.0, 0
        for _ in range(100):
            due, elapsed = F.periodic_update_due(F.DT_S, period, elapsed)
            fires += int(due)
            self.assertTrue(0.0 <= elapsed < period)
        # DT_S is the binary32 5 ms, marginally below 0.005, so 100 samples span
        # just under 0.5 s.  The scheduler is binary32, where the tolerance is
        # 16 * 2**-23 = 1.9e-6 rather than the 3.6e-15 of double, and that
        # tolerance absorbs the shortfall: the 100th sample fires, giving 25.
        # In a double-precision transcription the same word yields 24.  The
        # difference is not cosmetic -- it is the same tolerance that decides
        # S-firing starvation.
        self.assertEqual(fires, 25)
        due, _ = F.periodic_update_due(F.DT_S, 0.02, 0.0)
        self.assertFalse(due)

    def test_commit_period_rebases_the_timer(self):
        # set_pseudo_update_period_s applies elapsed = fmod(elapsed, new_period),
        # so a pending timer cannot survive a source commit unchanged.
        # Binary32: fmod(f32(0.13), f32(0.05)) is 0.0299999937, not 0.03.
        self.assertAlmostEqual(F.commit_period(0.13, 0.05), 0.03, places=7)
        self.assertEqual(F.commit_period(0.13, 0.05), F.f32(F.commit_period(0.13, 0.05)))
        self.assertEqual(F.commit_period(0.0, 0.05), 0.0)
        with self.assertRaises(RuntimeError):
            F.commit_period(float("nan"), 0.05)


class WordPropagation(unittest.TestCase):
    def test_zero_start_word_is_positive_and_covers_the_horizon(self):
        segs = [(1.1, 0.5, 10.0, F.WORD_SAMPLES)]
        P, fires, used = F.run_word(segs, [0.0, 0.0, 0.0, 0.0])
        self.assertEqual(used, F.WORD_SAMPLES)
        self.assertGreater(fires, 0)
        for i in range(4):
            self.assertTrue(math.isfinite(P[i][i]) and P[i][i] > 0.0)

    def test_zero_start_stays_below_a_seeded_start(self):
        # Riccati monotonicity in the initial covariance: P0 = 0 is a valid
        # lower below every admissible PSD initial covariance.
        segs = [(1.1, 0.5, 10.0, F.WORD_SAMPLES)]
        lo, _, _ = F.run_word(segs, [0.0, 0.0, 0.0, 0.0])
        hi, _, _ = F.run_word(segs, [1.0e6, 1.0e6, 1.0e6, 0.25])
        for i in range(4):
            self.assertLessEqual(lo[i][i], hi[i][i] * (1.0 + 1.0e-9))

    def test_a_changing_source_word_does_not_reuse_a_fixed_cell_firing_count(self):
        # The tau-scaled cadence and the fmod rebase mean a mixed word's firing
        # count is not any single cell's count.
        fixed_slow = F.run_word([(12.0, 0.05, 400.0, F.WORD_SAMPLES)], [0.0] * 4)[1]
        fixed_fast = F.run_word([(0.3333, 0.05, 400.0, F.WORD_SAMPLES)], [0.0] * 4)[1]
        mixed = F.run_word(
            [(12.0, 0.05, 400.0, 100), (0.3333, 0.05, 400.0, 100)] * 4, [0.0] * 4)[1]
        self.assertNotEqual(mixed, fixed_slow)
        self.assertNotEqual(mixed, fixed_fast)
        self.assertLess(fixed_slow, mixed)
        self.assertLess(mixed, fixed_fast)


class DiagnosticContract(unittest.TestCase):
    def test_build_is_non_promoting_and_validates(self):
        d = F.build(stride=211)
        self.assertEqual(validate_ok(d), [])
        self.assertTrue(d["non_promoting"])
        self.assertFalse(d["certifies_theorem_stage"])
        self.assertFalse(d["interval_certified"])
        self.assertFalse(d["quantifies_over_legal_histories"])
        self.assertEqual(d["canonical_gate"], 1.0e-18)
        self.assertEqual(d["word_samples"], 635)

    def test_validate_rejects_a_promoted_or_gate_moved_artifact(self):
        d = F.build(stride=211)
        for key, bad in (("non_promoting", False),
                         ("certifies_theorem_stage", True),
                         ("interval_certified", True),
                         ("canonical_gate", 1.0e-12)):
            broken = dict(d)
            broken[key] = bad
            self.assertTrue(F.validate(broken),
                            f"validate() accepted a tampered {key}")

    def test_every_row_forgets_its_initial_covariance(self):
        d = F.build(stride=211)
        for r in d["rows"]:
            self.assertTrue(r["P0_independent"],
                            f"node {r['source_node']} spread {r['P0_probe_relative_spread']}")


def validate_ok(d):
    return F.validate(d)


if __name__ == "__main__":
    unittest.main()


class AdversarialDescent(unittest.TestCase):
    """Single-cell dwell is not the adversary; the search must be able to say so."""

    def test_descent_finds_a_cheaper_legal_word_and_reports_convergence(self):
        # 0 -> {0,1}, 1 -> {0,1}, 2 unreachable. Node 1 scores lower.
        graph = {0: {0, 1}, 1: {0, 1}}
        score = {0: 1.0, 1: 0.5}

        def successors(n):
            return graph[n]

        def ratio_fn(seq):
            return min(score[n] for n in seq)

        out = F.adversarial_descent([0, 0, 0], successors, ratio_fn)
        self.assertTrue(out["converged"])
        self.assertLess(out["ratio"], 1.0)
        self.assertTrue(out["swaps"])
        # Every retained word must stay legal under the injected graph.
        for a, b in zip(out["word"], out["word"][1:]):
            self.assertIn(b, successors(a))

    def test_descent_respects_legality_and_can_find_nothing(self):
        # A graph with only self-loops admits no swap at all.
        graph = {0: {0}, 1: {1}}
        out = F.adversarial_descent([0, 0, 0], lambda n: graph[n], lambda s: 1.0)
        self.assertTrue(out["converged"])
        self.assertEqual(out["swaps"], [])
        self.assertEqual(out["word"], [0, 0, 0])

    def test_history_is_non_increasing(self):
        graph = {0: {0, 1, 2}, 1: {0, 1, 2}, 2: {0, 1, 2}}
        score = {0: 1.0, 1: 0.7, 2: 0.4}
        out = F.adversarial_descent([0, 0, 0], lambda n: graph[n],
                                    lambda s: min(score[n] for n in s))
        self.assertEqual(out["history"], sorted(out["history"], reverse=True))


class StarvationWitness(unittest.TestCase):
    """PR #476's zero-S word must be reproducible, or this tool is not conservative."""

    # Certified in #476's artifact ou3-p3-tau-ema-scheduler-cycle-diagnostic.
    TAU_LOW = 9.533334732055664
    TAU_HIGH = 9.533475875854492
    T_S_LOW = 0.1300000101327896
    T_S_HIGH = 0.13000193238258362

    def test_cadence_matches_the_certified_binary32_periods(self):
        # A double-precision 0.015/1.1 ratio misses these by ~1.5e-8 s, which is
        # enough to decide starvation against a 1.9e-6 tolerance.
        self.assertEqual(F.cadence_s(self.TAU_LOW), self.T_S_LOW)
        self.assertEqual(F.cadence_s(self.TAU_HIGH), self.T_S_HIGH)

    def test_alternating_tau_starves_S_while_holding_it_does_not(self):
        held = F.starvation_probe(self.TAU_LOW, self.TAU_LOW, 0.05, 400.0, gap=13)
        self.assertFalse(held["starved"])
        self.assertGreater(held["S_firings"], 0)
        alt = F.starvation_probe(self.TAU_LOW, self.TAU_HIGH, 0.05, 400.0, gap=13)
        self.assertTrue(alt["starved"])
        self.assertEqual(alt["S_firings"], 0)

    def test_a_starved_word_never_forgets_its_initial_covariance(self):
        alt = F.starvation_probe(self.TAU_LOW, self.TAU_HIGH, 0.05, 400.0, gap=13)
        self.assertFalse(alt["P0_independent"])
        self.assertGreater(alt["P0_probe_relative_spread"], 0.5)

    def test_build_records_that_the_dwell_result_is_conditional(self):
        d = F.build(stride=401)
        self.assertTrue(d["dwell_result_conditional_on_S_firing"])
