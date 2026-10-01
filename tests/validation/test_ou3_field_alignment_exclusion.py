"""Analytical substitution regressions; not trajectory or IEEE-754 certification."""
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))
from ou3_theorem.field_alignment_exclusion import (  # noqa: E402
    geometric_margin, physical_mean_charge, sin_ten_degrees_interval,
    small_branch_process_bound, trapezoidal_weights,
)


class FieldAlignmentExclusionTests(unittest.TestCase):
    def test_literal_polynomial_uniform_bound(self):
        bound = small_branch_process_bound()
        self.assertLess(bound, F(1024, 1000))
        self.assertLess(bound, F(103, 100))

    def test_safe_covariance_radius_coefficient(self):
        self.assertGreater(F(203, 50)**2, 16 * F(103, 100))

    def test_trigonometric_enclosure(self):
        lo, hi = sin_ten_degrees_interval()
        self.assertGreater(lo, F("0.1736481776669303488"))
        self.assertLess(hi, F("0.1736481776669303490"))
        self.assertLess(lo, hi)

    def test_nonuniform_weights_and_effective_gap(self):
        gaps = (F(1, 250), F(1, 200), F(3, 500))
        weights = trapezoidal_weights(gaps)
        self.assertEqual(sum(weights), 1)
        self.assertEqual(len(weights), len(gaps) + 1)
        effective = sum(gap**2 for gap in gaps) / sum(gaps)
        self.assertLessEqual(effective, max(gaps))

    def test_reject_nonpositive_gaps(self):
        for gaps in ((), (F(0),), (F(-1), F(2))):
            with self.assertRaises(ValueError):
                trapezoidal_weights(gaps)

    def test_correct_accelerometer_cadence_charge(self):
        charge = physical_mean_charge(duration=F(100), velocity_bound=F(11, 2),
                                      jerk_bound=F(100), effective_gap=F(3, 500))
        self.assertEqual(charge, F(26, 100))
        old = physical_mean_charge(duration=F(100), velocity_bound=F(11, 2),
                                   jerk_bound=F(100), effective_gap=F(1, 25))
        self.assertEqual(old - charge, F(85, 100))

    def margin(self, radius, defect):
        lo, _ = sin_ten_degrees_interval()
        return geometric_margin(
            transverse_gravity_lower=F("9.80665") * lo, duration=F(100),
            velocity_bound=F(11, 2), jerk_bound=F(100), effective_gap=F(3, 500),
            acceleration_error_factor=F(203, 50), retained_radius=radius,
            force_defect=defect,
        )

    def test_positive_conditional_margin(self):
        self.assertGreater(self.margin(F(1, 4), F(0)), F("0.4279"))

    def test_actual_duration_after_word_boundary_trimming(self):
        lo, _ = sin_ten_degrees_interval()
        value = geometric_margin(
            transverse_gravity_lower=F("9.80665") * lo, duration=F("99.988"),
            velocity_bound=F(11, 2), jerk_bound=F(100), effective_gap=F(3, 500),
            acceleration_error_factor=F(203, 50), retained_radius=F(1, 4),
            force_defect=F(0),
        )
        self.assertGreater(value, F("0.4278"))
        self.assertLess(self.margin(F(1, 4), F(0)) - value, F("0.000014"))

    def test_defects_are_not_silently_discarded(self):
        self.assertLess(self.margin(F(1, 4), F(1, 2)), 0)

    def test_larger_radius_not_promoted(self):
        self.assertLess(self.margin(F(1), F(0)), 0)

    def test_startup_fraction_only_fixed_reference_specialization(self):
        margin = geometric_margin(
            transverse_gravity_lower=F("9.80665") / 20, duration=F(100),
            velocity_bound=F(11, 2), jerk_bound=F(100), effective_gap=F(3, 500),
            acceleration_error_factor=F(203, 50), retained_radius=F(1, 25),
            force_defect=F(0),
        )
        self.assertGreater(margin, F("0.0679"))

    def test_sharp_scalar_sampling_cell(self):
        # a(s)=c-J*min(s,h-s) saturates the endpoint integration inequality.
        c, jerk, gap = F(2), F(100), F(3, 500)
        integral = c*gap - jerk*gap**2/4
        trapezoidal = c*gap
        self.assertEqual(trapezoidal - integral, jerk*gap**2/4)

    def test_negative_defect_is_rejected(self):
        with self.assertRaises(ValueError):
            self.margin(F(1, 4), F(-1))


if __name__ == "__main__":
    unittest.main()
