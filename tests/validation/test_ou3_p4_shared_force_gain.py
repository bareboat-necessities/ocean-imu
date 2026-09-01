from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p4_candidate_first_accel_range_v3 as V3
import ou3_p4_shared_force_gain as SHARED


def _brute_force_sup(fn, m: Interval, samples: int = 20001) -> float:
    lo, hi = m.lo, m.hi
    best = 0.0
    for i in range(samples):
        mv = lo + (hi - lo) * i / (samples - 1)
        best = max(best, fn(mv))
    return best


class SharedForceGainMathTests(unittest.TestCase):
    def test_self_p_row_bounds_the_true_supremum(self):
        for m in (Interval(5.0, 7.83), Interval(12.2, 19.2), Interval(0.9, 40.0)):
            for p in (1.225e-3, 4.0e-3, 0.5):
                for lam in (3.3e-3, 0.29, 36.0):
                    bound = SHARED.sup_self_p_gain(m, p, lam)
                    truth = _brute_force_sup(lambda mv: mv * p / (mv * mv * p + lam), m)
                    self.assertGreaterEqual(bound, truth)
                    self.assertLess(bound, truth * 1.0001 + 1.0e-12)

    def test_independent_row_bounds_the_true_supremum(self):
        for m in (Interval(5.0, 7.83), Interval(19.2, 30.0)):
            for p in (1.225e-3, 4.0e-3):
                for lam in (3.3e-3, 0.29):
                    for c in (1.5e-3, 2.1e-3):
                        bound = SHARED.sup_independent_gain(m, c, p, lam)
                        truth = _brute_force_sup(
                            lambda mv: mv * c / (mv * mv * p + lam), m)
                        self.assertGreaterEqual(bound, truth)
                        self.assertLess(bound, truth * 1.0001 + 1.0e-12)

    def test_interior_maximum_is_used_when_reachable(self):
        # m* = sqrt(lam/p) = 10 lies inside the cell, so the endpoint values
        # alone would understate the supremum.
        p, lam = 1.0e-2, 1.0
        m = Interval(2.0, 50.0)
        bound = SHARED.sup_self_p_gain(m, p, lam)
        self.assertGreaterEqual(bound, 0.5 * math.sqrt(p / lam))
        truth = _brute_force_sup(lambda mv: mv * p / (mv * mv * p + lam), m)
        self.assertGreaterEqual(bound, truth)

    def test_degenerate_domains_are_rejected(self):
        with self.assertRaises(RuntimeError):
            SHARED.sup_self_p_gain(Interval(0.0, 1.0), 1.0e-3, 1.0e-3)
        with self.assertRaises(RuntimeError):
            SHARED.sup_self_p_gain(Interval(1.0, 2.0), 1.0e-3, 0.0)
        with self.assertRaises(RuntimeError):
            SHARED.sup_independent_gain(Interval(1.0, 2.0), 1.0e-3, 0.0, 1.0e-3)


class SharedForceGainRowTests(unittest.TestCase):
    CELL = dict(
        tilt=0.0012250000000000002,
        yaw=0.007568999999999999,
        eps=1.4828249737599615e-10,
        x=Interval(0.0, 0.0625),
        m=Interval(5.0, 7.825422900366438),
        paw=Interval(0.0024261138354981, 0.30093178099430606),
        racc_var=Interval(0.0008666185998411353, 0.0008666185998411355),
    )

    def test_gain_rows_are_tighter_than_the_naive_interval_rows(self):
        k_naive, kh_naive, dn = V3._tangent_structured_gain_bounds(**self.CELL)
        k_exact, kh_exact, de = SHARED.shared_force_structured_gain_bounds(**self.CELL)
        self.assertLess(k_exact, k_naive)
        self.assertLessEqual(kh_exact, kh_naive)
        for row in ("g_perp_upper", "g_u_upper", "g_z_upper"):
            self.assertLessEqual(de[row], dn[row])
        self.assertTrue(de["shared_force_magnitude_dependency_preserved"])
        self.assertGreater(k_naive / k_exact, 1.5)

    def test_psd_remainder_and_KH_treatment_are_unchanged(self):
        _k, _kh, dn = V3._tangent_structured_gain_bounds(**self.CELL)
        _ke, _khe, de = SHARED.shared_force_structured_gain_bounds(**self.CELL)
        for key in ("attitude_PSD_remainder_operator_upper",
                    "nominal_tangent_innovation_lower",
                    "tangent_innovation_perturbation_upper",
                    "perturbed_tangent_innovation_lower",
                    "PSD_innovation_perturbation_tangent_only",
                    "PSD_innovation_axial_row_column_exact_zero"):
            self.assertEqual(de[key], dn[key])


class SharedForceGainCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = SHARED.build()

    def test_lemma_validates_and_is_uniformly_tighter(self):
        self.assertEqual(SHARED.validate(self.d), [])
        self.assertTrue(self.d["uniformly_at_least_as_tight_as_interval_gain"])
        self.assertGreaterEqual(self.d["minimum_observed_tightening_factor"], 1.0)
        # The amount of improvement varies with the source covariance (notably
        # horizontal R_S). Require genuine improvement, not a tuning-dependent
        # 1.5x diagnostic threshold; all-cell non-regression is checked above.
        self.assertGreater(self.d["maximum_observed_tightening_factor"], 1.0)

    def test_lemma_promotes_nothing(self):
        self.assertFalse(self.d["filter_changed"])
        self.assertFalse(self.d["source_replay_used"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(self.d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])


if __name__ == "__main__":
    unittest.main()
