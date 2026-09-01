from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p4_first_accel_aw_sigma_consistency as CONSIST


class ForceSubdivisionTests(unittest.TestCase):
    def test_subcells_tile_the_cell_and_keep_its_endpoints(self):
        m = Interval(5.0, 30.0)
        for pieces in (1, 2, 8, 16):
            cells = CONSIST._force_subcells(m, pieces)
            self.assertEqual(len(cells), pieces)
            self.assertLessEqual(cells[0].lo, m.lo)
            self.assertGreaterEqual(cells[-1].hi, m.hi)
            for a, b in zip(cells, cells[1:]):
                self.assertGreaterEqual(b.hi, a.hi)
                self.assertLessEqual(a.lo, b.lo)

    def test_a_nonpositive_subdivision_is_rejected(self):
        with self.assertRaises(ValueError):
            CONSIST._force_subcells(Interval(5.0, 30.0), 0)


class ConsistencySearchTests(unittest.TestCase):
    """The bisection is exercised on a synthetic table with a known answer."""

    def _table(self, k, sigma_hi, aw, ba):
        return {
            "rows": [{
                "aw_after_prefix_upper_mps2": aw,
                "tuner_sigma_aw_upper_mps2": sigma_hi,
                "force_upper_mps2": 5.0,
                "Ktheta_norm_upper": k,
            }],
            "accel_bias_error_norm_upper_mps2": ba,
        }

    def test_worst_nuisance_clamps_the_aw_error_at_c_times_sigma(self):
        t = self._table(k=0.3, sigma_hi=0.5, aw=3.0, ba=0.5)
        loose, _w = CONSIST._worst_nuisance(t, 0.0, math.inf)
        tight, _w = CONSIST._worst_nuisance(t, 0.0, 2.0)
        self.assertAlmostEqual(loose, 0.3 * (3.0 + 0.5), delta=1e-9)
        self.assertAlmostEqual(tight, 0.3 * (1.0 + 0.5), delta=1e-9)
        self.assertLess(tight, loose)

    def test_a_larger_c_never_lowers_the_nuisance(self):
        t = self._table(k=0.3, sigma_hi=0.5, aw=3.0, ba=0.5)
        vals = [CONSIST._worst_nuisance(t, 0.0, c)[0] for c in (0.0, 0.5, 1.0, 4.0, 12.0)]
        self.assertEqual(vals, sorted(vals))

    def test_critical_constant_brackets_the_budget_crossing(self):
        t = self._table(k=0.3, sigma_hi=0.5, aw=3.0, ba=0.5)
        # zero-aw residual is 0.15; budget 0.30 leaves room for 0.5 m/s^2 of
        # a_w error, i.e. c = 1.0.
        out = CONSIST._critical_consistency_constant(t, 0.0, 0.30)
        self.assertTrue(out["any_finite_constant_closes_this_angle"])
        self.assertAlmostEqual(out["critical_consistency_constant"], 1.0, delta=1e-6)
        lo, hi = out["bracket"]
        self.assertLess(CONSIST._worst_nuisance(t, 0.0, lo)[0], 0.30)
        self.assertGreaterEqual(CONSIST._worst_nuisance(t, 0.0, hi)[0], 0.30)

    def test_no_constant_helps_when_the_zero_error_residual_already_exceeds(self):
        t = self._table(k=0.3, sigma_hi=0.5, aw=3.0, ba=2.0)
        out = CONSIST._critical_consistency_constant(t, 0.0, 0.30)
        self.assertFalse(out["any_finite_constant_closes_this_angle"])
        self.assertEqual(out["critical_consistency_constant"], 0.0)
        self.assertGreaterEqual(out["nuisance_at_zero_aw_error_rad"], 0.30)

    def test_an_already_fitting_pairing_reports_an_unbounded_constant(self):
        t = self._table(k=0.01, sigma_hi=0.5, aw=3.0, ba=0.5)
        out = CONSIST._critical_consistency_constant(t, 0.0, 0.30)
        self.assertTrue(out["unconstrained_pairing_already_fits"])
        self.assertEqual(out["critical_consistency_constant"], math.inf)

    def test_a_nonpositive_budget_admits_no_constant(self):
        t = self._table(k=0.3, sigma_hi=0.5, aw=3.0, ba=0.5)
        out = CONSIST._critical_consistency_constant(t, 0.0, 0.0)
        self.assertFalse(out["any_finite_constant_closes_this_angle"])


class ConsistencyCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CONSIST.build()

    def test_certificate_validates_and_promotes_nothing(self):
        d = self.d
        self.assertEqual(CONSIST.validate(d), [])
        self.assertTrue(d["distance_only_no_verdict_emitted"])
        self.assertTrue(d["consistency_constant_is_a_conditional_requirement_not_a_theorem"])
        self.assertFalse(d["aw_sigma_consistency_declared_in_domain"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["declared_entrance_shrunk"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])

    def test_joint_force_pairing_is_an_unconditional_tightening(self):
        self.assertTrue(self.d["joint_force_pairing_is_unconditional"])
        for r in self.d["angle_rows"]:
            self.assertLessEqual(r["nuisance_joint_force_pairing_upper_rad"],
                                 r["nuisance_separate_force_bounds_upper_rad"])
            self.assertGreaterEqual(r["joint_force_pairing_tightening_factor"], 1.0)

    def test_the_unconstrained_pairing_still_fits_nowhere(self):
        for r in self.d["angle_rows"]:
            self.assertFalse(r["joint_pairing_fits_inside_budget"])
            self.assertGreater(r["nuisance_over_budget_ratio_joint"], 1.0)

    def test_the_zero_aw_residual_alone_closes_the_wide_candidates_out(self):
        rows = {r["angle_deg"]: r for r in self.d["angle_rows"]}
        for deg in (30.0, 35.0, 40.0, 45.0):
            r = rows[deg]
            self.assertGreater(r["nuisance_at_zero_aw_error_rad"],
                               r["sector_invariance_correction_budget_upper_rad"])
            self.assertFalse(r["any_finite_constant_closes_this_angle"])
            self.assertEqual(r["critical_consistency_constant"], 0.0)

    def test_the_narrow_candidates_close_only_with_a_finite_constant(self):
        rows = {r["angle_deg"]: r for r in self.d["angle_rows"]}
        for deg in (15.0, 20.0, 25.0):
            r = rows[deg]
            self.assertTrue(r["any_finite_constant_closes_this_angle"])
            self.assertGreater(r["critical_consistency_constant"], 0.0)
            self.assertLess(r["critical_consistency_constant"],
                            r["unconstrained_c_at_worst_cell"])
            self.assertGreater(r["consistency_tightening_needed_factor"], 1.0)
        self.assertEqual(self.d["widest_angle_closed_by_a_finite_consistency_constant"], 25.0)

    def test_the_required_constant_relaxes_as_the_candidate_narrows(self):
        rows = sorted((r for r in self.d["angle_rows"]
                       if r["any_finite_constant_closes_this_angle"]),
                      key=lambda r: r["angle_deg"])
        cs = [r["critical_consistency_constant"] for r in rows]
        self.assertEqual(cs, sorted(cs, reverse=True))

    def test_the_tuner_sigma_range_is_the_deployed_safety_envelope(self):
        lo, hi = self.d["tuner_sigma_aw_range_mps2"]
        self.assertAlmostEqual(lo, 0.05, delta=1e-9)
        self.assertAlmostEqual(hi, 6.0, delta=1e-9)


if __name__ == "__main__":
    unittest.main()
