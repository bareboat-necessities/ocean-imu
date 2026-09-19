"""The H18 service-superword feasibility diagnostic measures; it never promotes.

These tests fix three things: the high-precision linear algebra returns what it
claims, a degenerate or non-H18 export is refused rather than evaluated, and no
outcome of the diagnostic -- including a synthetic contracting one -- discharges
a proof obligation.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem import h18_superword as diag

N = diag.ERROR_DIMENSION
STEP = -0.002   # the exported root response is a measured displacement, not +1


def rows(matrix):
    return [value for row in matrix for value in row]


def scaled_identity(scale):
    return [[scale if i == j else 0.0 for j in range(N)] for i in range(N)]


def diagonal(values):
    return [[values[i] if i == j else 0.0 for j in range(N)] for i in range(N)]


def product(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(N)) for j in range(N)]
            for i in range(N)]


def upper(matrix):
    return [matrix[i][j] for i in range(N) for j in range(i, N)]


def synthetic_export(gain=0.5, samples=4, covariance=1.0, mag_prefix=2):
    """An export whose superword map is a known power of a diagonal gain."""
    transition = diagonal([gain] * N)
    difference = [[STEP if i == j else 0.0 for j in range(N)] for i in range(N)]
    differences = [difference]
    for _ in range(samples):
        difference = product(transition, difference)
        differences.append(difference)
    metric = upper(scaled_identity(covariance))
    return {
        "error_dimension": N,
        "superword_samples": samples,
        "root_time_s": 12.5,
        "reference_error": [[0.01] * N for _ in range(samples + 1)],
        "covariance_upper": [list(metric) for _ in range(samples + 1)],
        "error_difference": [rows(m) for m in differences],
        "root_difference_coarse": rows(differences[0]),
        "endpoint_difference_coarse": rows(differences[-1]),
        "root_rotation_world_to_body": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        "applied_magnetic_corrections": [{
            "prefix": mag_prefix,
            "time_s": 12.51,
            "sensitivity_axis": [20.0, 0.0, 43.0],
            "innovation_covariance": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            "phase": "pre_correction",
            "pre_error_difference": rows(differences[0]),
        }],
        "acc_bias_held_through_superword": True,
        "attitude_injection_finite": True,
        "inherited_state": True,
    }


class LinearAlgebraTests(unittest.TestCase):
    def test_storage_matches_the_inverse_metric(self):
        with localcontext() as context:
            context.prec = 50
            covariance = [[Decimal(4), Decimal(1)], [Decimal(1), Decimal(3)]]
            factor = diag.ldl_factor(covariance)
            e = [Decimal(2), Decimal(-1)]
            # P^-1 = [[3,-1],[-1,4]]/11, so e^T P^-1 e = (12 + 4 + 4)/11.
            self.assertAlmostEqual(float(diag.storage(factor, e)), 20.0 / 11.0, places=12)

    def test_indefinite_metric_is_refused(self):
        with self.assertRaises(diag.SuperwordExportError):
            diag.ldl_factor([[Decimal(1), Decimal(2)], [Decimal(2), Decimal(1)]])

    def test_inverse_is_exact_enough_to_round_trip(self):
        with localcontext() as context:
            context.prec = 50
            a = [[Decimal(2), Decimal(1), Decimal(0)],
                 [Decimal(1), Decimal(3), Decimal(1)],
                 [Decimal(0), Decimal(1), Decimal(4)]]
            inv = diag.inverse(a)
            for i in range(3):
                for j in range(3):
                    entry = sum(a[i][k] * inv[k][j] for k in range(3))
                    self.assertAlmostEqual(float(entry), 1.0 if i == j else 0.0, places=30)

    def test_singular_root_response_is_refused(self):
        with self.assertRaises(diag.SuperwordExportError):
            diag.inverse([[Decimal(1), Decimal(2)], [Decimal(2), Decimal(4)]])


class WorstRatioTests(unittest.TestCase):
    def test_known_diagonal_map_gives_its_squared_gain(self):
        report = diag.evaluate(synthetic_export(gain=0.5, samples=4), precision=40)
        # Psi = 0.5^4 I in an identity metric, so the worst ratio is (0.5^4)^2.
        self.assertAlmostEqual(report["worst_admissible_ratio"], 0.5 ** 8, places=12)
        self.assertTrue(report["strict_contraction_observed"])
        self.assertTrue(report["endpoint_direction_prefix_nonexpansive"])

    def test_the_worst_direction_is_the_largest_gain(self):
        export = synthetic_export(gain=0.5, samples=1, mag_prefix=1)
        gains = [0.5] * N
        gains[9] = 0.9
        difference = [[STEP if i == j else 0.0 for j in range(N)] for i in range(N)]
        export["error_difference"] = [
            rows(difference), rows(product(diagonal(gains), difference))]
        export["endpoint_difference_coarse"] = export["error_difference"][-1]
        report = diag.evaluate(export, precision=40)
        self.assertAlmostEqual(report["worst_admissible_ratio"], 0.81, places=12)
        blocks = report["limiting_direction_blocks"]
        self.assertGreater(blocks["position"], 0.999)

    def test_a_dominant_direction_no_seed_favours_is_still_found(self):
        # The only nonzero gain sits on the last coordinate, so a map that an
        # unlucky single seed would miss entirely. The seeds span the space, so
        # the reported ratio is the real maximum rather than a subdominant one.
        export = synthetic_export(gain=0.0, samples=1, mag_prefix=1)
        gains = [0.0] * N
        gains[N - 1] = 2.0
        difference = [[STEP if i == j else 0.0 for j in range(N)] for i in range(N)]
        export["error_difference"] = [
            rows(difference), rows(product(diagonal(gains), difference))]
        export["endpoint_difference_coarse"] = export["error_difference"][-1]
        report = diag.evaluate(export, precision=40)
        self.assertAlmostEqual(report["worst_admissible_ratio"], 4.0, places=12)
        self.assertFalse(report["strict_contraction_observed"])
        self.assertGreater(report["limiting_direction_blocks"]["accelerometer_bias"], 0.999)

    def test_an_unconverged_ratio_fails_closed(self):
        export = synthetic_export(gain=0.5, samples=1, mag_prefix=1)
        with localcontext() as context:
            context.prec = 40
            n = N
            differences = [diag.square_from_rows(row, n)
                           for row in export["error_difference"]]
            covariances = [diag.symmetric_from_upper(row, n)
                           for row in export["covariance_upper"]]
            superword = diag.SuperwordMap(differences[-1], diag.inverse(differences[0]))
            metric = diag.ldl_factor(covariances[0])
            with self.assertRaises(diag.SuperwordExportError):
                diag.worst_admissible_ratio(superword, covariances[0], metric, metric,
                                            n, iterations=1)

    def test_the_reported_eigenpair_carries_its_residual(self):
        report = diag.evaluate(synthetic_export(gain=0.5, samples=2), precision=40)
        self.assertLess(report["eigenpair_relative_residual"], 1e-20)
        self.assertEqual(report["power_iteration_seeds"], N)

    def test_growth_is_reported_without_being_promoted(self):
        report = diag.evaluate(synthetic_export(gain=1.4, samples=2), precision=40)
        self.assertGreater(report["worst_admissible_ratio"], 1.0)
        self.assertFalse(report["strict_contraction_observed"])
        self.assertFalse(report["endpoint_direction_prefix_nonexpansive"])
        self.assertFalse(report["certificate_complete"])


class MagneticInformationTests(unittest.TestCase):
    def test_untransported_service_is_rank_deficient(self):
        # With an identity map the axial gyro-bias coordinate never reaches the
        # magnetometer, so frequent valid corrections still leave the pair
        # unobservable. A maximum event gap alone cannot supply service.
        report = diag.evaluate(synthetic_export(gain=1.0, samples=2), precision=40)
        self.assertAlmostEqual(report["magnetic_information_min_eigenvalue"], 0.0, places=20)
        self.assertFalse(report["magnetic_information_meets_floor"])

    def test_applied_corrections_are_counted_from_the_export_only(self):
        report = diag.evaluate(synthetic_export(), precision=40)
        self.assertEqual(report["applied_magnetic_corrections"], 1)


class ValidationTests(unittest.TestCase):
    def test_required_fields_are_enforced(self):
        for key in diag.REQUIRED_KEYS:
            export = synthetic_export()
            del export[key]
            with self.subTest(field=key), self.assertRaises(diag.SuperwordExportError):
                diag.validate(export)

    def test_a_superword_without_an_applied_correction_is_refused(self):
        export = synthetic_export()
        export["applied_magnetic_corrections"] = []
        with self.assertRaises(diag.SuperwordExportError):
            diag.validate(export)

    def test_a_released_or_uninherited_run_is_not_an_h18_superword(self):
        for flag in ("acc_bias_held_through_superword", "attitude_injection_finite",
                     "inherited_state"):
            export = synthetic_export()
            export[flag] = False
            with self.subTest(flag=flag), self.assertRaises(diag.SuperwordExportError):
                diag.validate(export)

    def test_a_correction_outside_the_window_cannot_inform_it(self):
        for prefix in (0, 99):
            export = synthetic_export(mag_prefix=prefix)
            with self.subTest(prefix=prefix), self.assertRaises(diag.SuperwordExportError):
                diag.validate(export)

    def test_truncated_prefixes_cannot_pass(self):
        export = synthetic_export()
        export["error_difference"] = export["error_difference"][:-1]
        with self.assertRaises(diag.SuperwordExportError):
            diag.validate(export)

    def test_empty_and_nonfinite_exports_are_refused(self):
        export = synthetic_export()
        export["superword_samples"] = 0
        with self.assertRaises(diag.SuperwordExportError):
            diag.validate(export)
        export = synthetic_export()
        export["reference_error"][0][0] = float("inf")
        with self.assertRaises(diag.SuperwordExportError):
            diag.validate(export)


class NonPromotionTests(unittest.TestCase):
    def test_no_outcome_discharges_an_obligation(self):
        for gain in (0.2, 1.0, 1.9):
            report = diag.evaluate(synthetic_export(gain=gain), precision=40)
            with self.subTest(gain=gain):
                self.assertFalse(report["certificate_complete"])
                self.assertFalse(report["obligation_discharged"])
                self.assertFalse(report["source_uniform"])
                self.assertFalse(report["supply_constant_evaluated"])
                self.assertTrue(report["local_incremental_map_only"])

    def test_the_held_bias_obstruction_travels_with_the_report(self):
        report = diag.evaluate(synthetic_export(gain=1.0), precision=40)
        obstruction = report["held_bias_obstruction"]
        # A unit-gain synthetic map reproduces the held bias and keeps the
        # covariance block fixed, which is exactly the obstruction's premise.
        self.assertTrue(obstruction["non_contraction_obstruction"])
        self.assertEqual(obstruction["full_state_rho_status"], "excluded")
        self.assertFalse(obstruction["obstruction_is_an_instability_claim"])

    def test_a_moving_bias_block_is_not_claimed_as_an_obstruction(self):
        report = diag.evaluate(synthetic_export(gain=0.5), precision=40)
        obstruction = report["held_bias_obstruction"]
        self.assertFalse(obstruction["non_contraction_obstruction"])
        self.assertEqual(obstruction["full_state_rho_status"], "undecided_here")


if __name__ == "__main__":
    unittest.main()
