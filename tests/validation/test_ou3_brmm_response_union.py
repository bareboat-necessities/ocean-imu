from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/stability"))
import ou3_brmm_response_union as U  # noqa: E402


class ResponseUnionTests(unittest.TestCase):
    def test_union_keeps_linear_and_stokes_distinct(self):
        d = U.build()
        self.assertEqual(U.validate(d), [])
        self.assertTrue(d["linear_vessel_domain_unchanged"])
        self.assertFalse(d["legacy_continuum_replaced_by_finite_grid"])
        self.assertEqual(d["linear_RAO_moment_and_surface_period_lemmas_scope"], [U.LINEAR])

    def test_per_event_branch_switch_or_missing_branch_rejected(self):
        for key, value in (("per_event_branch_switching_allowed", True),
                           ("branches", [U.STOKES])):
            d = deepcopy(U.build())
            d[key] = value
            self.assertTrue(U.validate(d))

    def test_exact_stokes_bounds_retain_higher_harmonics(self):
        self.assertEqual(U.stokes_hard_bounds([("2", "1", "1/10")]), {
            "B0_exact": "843/8000", "B1_exact": "889/4000",
            "B2_exact": "987/2000", "slope_bound_exact": "889/8000"})

    def test_steepness_is_model_premise_not_independent_coefficient(self):
        U.stokes_hard_bounds([("2", "1", "1/5")])
        with self.assertRaises(ValueError):
            U.stokes_hard_bounds([("2", "1", "201/1000")])

    def test_partial_branch_or_mode_cannot_certify_union(self):
        d = U.branch_mode_coverage({"H18": True, "A21": True})
        self.assertTrue(U.all_branches_modes_closed(d))
        d[U.STOKES]["A21"] = False
        self.assertFalse(U.all_branches_modes_closed(d))
        del d[U.STOKES]
        self.assertFalse(U.all_branches_modes_closed(d))
        with self.assertRaises(ValueError):
            U.branch_mode_coverage({"H18": True})


if __name__ == "__main__":
    unittest.main()
