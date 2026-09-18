import sys
from dataclasses import replace
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem.finite_error import (
    ServiceSuperwordWitness, audit_service_superword, incomplete_diagnostic,
)


class FiniteErrorTargetTests(unittest.TestCase):
    def witness(self):
        return ServiceSuperwordWitness(
            history_id="h", V_start=1.0, V_end=0.8,
            disturbance_energy=0.0, rho_candidate=0.9, disturbance_gain=0.0,
            prefix_values=(1.0, 0.9, 0.8), prefix_bounds=(1.1, 1.0, 0.9),
            magnetic_information_min_eigenvalue=1.2, magnetic_information_floor=1.0,
            complete_shipping_map=True, inherited_state=True,
            same_history_motion=True, same_history_bias=True,
            applied_magnetic_information=True, finite_error_map=True,
            arithmetic_enclosed=True,
        )

    def test_point_diagnostic_remains_incomplete(self):
        w = incomplete_diagnostic("h", (1.0, 0.95, 0.90), 1.0)
        r = audit_service_superword(w)
        self.assertTrue(r["strict_rho"])
        self.assertIsNone(w.prefix_bounds)
        self.assertFalse(r["prefix_retention"])
        self.assertFalse(r["structural_premises"])
        self.assertFalse(r["magnetic_information"])
        self.assertFalse(r["certificate_complete"])

    def test_consistent_numbers_and_flags_are_not_a_certificate(self):
        r = audit_service_superword(self.witness())
        self.assertTrue(r["point_consistency_pass"])
        self.assertFalse(r["structural_premises_verified_here"])
        self.assertFalse(r["certificate_complete"])
        self.assertFalse(r["point_witness_promoted_to_source_uniform_theorem"])

    def test_missing_prefixes_cannot_pass_vacuously(self):
        with self.assertRaises(ValueError):
            replace(self.witness(), prefix_values=(), prefix_bounds=())
        self.assertFalse(audit_service_superword(
            replace(self.witness(), prefix_bounds=None))["point_consistency_pass"])

    def test_prefix_endpoints_must_match(self):
        with self.assertRaises(ValueError):
            replace(self.witness(), prefix_values=(0.0, 0.0, 0.0))

    def test_infinite_or_nan_bounds_are_rejected(self):
        for bad in (float("inf"), float("nan")):
            with self.subTest(value=bad), self.assertRaises(ValueError):
                replace(self.witness(), prefix_bounds=(1.1, bad, 0.9))

    def test_rhs_overflow_cannot_pass(self):
        w = replace(self.witness(), disturbance_energy=1e308, disturbance_gain=1e308)
        r = audit_service_superword(w)
        self.assertFalse(r["arithmetic_finite"])
        self.assertFalse(r["dissipativity"])
        self.assertFalse(r["point_consistency_pass"])
        self.assertIsNone(r["dissipativity_rhs"])

    def test_assertions_must_be_booleans(self):
        with self.assertRaises(ValueError):
            replace(self.witness(), arithmetic_enclosed="yes")

    def test_prefix_information_and_rho_failures_are_visible(self):
        for w in (
            replace(self.witness(), prefix_bounds=(0.1, 0.1, 0.1)),
            replace(self.witness(), magnetic_information_min_eigenvalue=0.0),
            replace(self.witness(), rho_candidate=1.0),
        ):
            self.assertFalse(audit_service_superword(w)["point_consistency_pass"])


if __name__ == "__main__":
    unittest.main()
