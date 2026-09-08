from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/stability"))
from ou3_interval import Interval, matrix_point, matrix_identity  # noqa: E402
import ou3_p4_bounded_bias_motion as M  # noqa: E402
import ou3_p4_complete_brmm_joint_sector_master as J  # noqa: E402


def contract():
    return M.build(p3_contract={
        "mems_bias_preconditions": M.BIAS.build(), "P3_CONDITIONAL_BRMM_PASS": True,
        "modes": {mode: {"relative_Riccati_injection_margin_lower": 1e-18}
                  for mode in ("H18", "A21")}})


class BoundedBiasMotionTests(unittest.TestCase):
    def test_clamp_bounds_estimate_not_error_or_true_bias(self):
        self.assertEqual(M.bias_error_bound("0.1"), F(1, 2))
        self.assertEqual(M.bias_error_bound("0.4"), F(4, 5))
        self.assertEqual(M.bias_error_bound("0.7"), F(11, 10))
        for bound in (-1, None, float("inf")):
            with self.assertRaises((ValueError, TypeError)):
                M.bias_error_bound(bound)

    def bounds(self, **overrides):
        args = dict(rho="1/2", Gamma=2, true_bias_bound="1/10", duration=3,
                    gamma=("1/10", 2, 1), kappa=("1/5", 1, 1),
                    source_energy="1/50", noise_energy=0, endpoint_level=1, chart_level=3)
        return M.compose_bounds(**(args | overrides))

    def test_geometric_series_and_transient_gain_are_exact(self):
        d = self.bounds()
        self.assertEqual(d["bias_window_energy_bound"], F(3, 4))
        self.assertEqual(d["endpoint_ultimate_bound"], F(23, 100))
        self.assertEqual(d["all_prefix_ultimate_bound"], F(63, 100))
        self.assertTrue(d["endpoint_level_invariant"])
        self.assertTrue(d["strict_prefix_chart_retention"])
        self.assertFalse(d["actual_motion_gains_certified_by_this_arithmetic"])
        w = F(1)
        for _ in range(20):
            w = w/2+d["endpoint_forcing_bound"]
        self.assertEqual(w, F(1, 2)**20+(1-F(1, 2)**20)*F(23, 100))

    def test_noise_off_does_not_remove_model_mismatch_error_floor(self):
        d = self.bounds(gamma=(0, 2, 1), kappa=(0, 1, 1), noise_energy=0)
        self.assertEqual(d["endpoint_ultimate_bound"], F(2, 25))
        self.assertEqual(d["all_prefix_ultimate_bound"], F(9, 50))
        self.assertGreater(d["all_prefix_ultimate_bound"], 0)

    def test_bounded_bias_does_not_automatically_retain_chart(self):
        self.assertFalse(self.bounds(chart_level=2)["strict_prefix_chart_retention"])
        self.assertFalse(self.bounds(endpoint_level="1/10")["endpoint_level_invariant"])
        for overrides in ({"rho": 1}, {"rho": 0}, {"Gamma": "1/2"},
                          {"gamma": (-1, 0, 0)}, {"duration": 0}):
            with self.assertRaises(ValueError):
                self.bounds(**overrides)

    def maps(self):
        c0 = [[float(i == j) for j in range(21)] for i in range(18)]
        cn = [[0.5*x for x in row] for row in c0]
        cn[0][18] = 1.0
        qb = [[float(i == j and i >= 18) for j in range(21)] for i in range(21)]
        return matrix_point(c0), matrix_point(cn), matrix_point(qb)

    def test_full_graph_master_keeps_bias_to_motion_cross_term(self):
        c0, cn, qb = self.maps()
        master = M.motion_master(c0, cn, matrix_identity(18), matrix_identity(18),
                                 Interval.point(.5), [(Interval.point(4), qb)])
        self.assertEqual(len(master), 21)
        self.assertTrue(master[0][18].contains(.5))
        self.assertTrue(master[18][18].contains(-3))
        z = [Interval.point(float(i in (0, 18))) for i in range(21)]
        self.assertTrue(J.quadratic_value(master, z).contains(1.5**2-.5-4))
        # Synthetic algebra witness only, never an Ocean-IMU gain certificate.
        zero_sector = matrix_point([[0.]*21 for _ in range(21)])
        passed, pivots = J.certify_strict_joint_sector_domination(master, [zero_sector], [0.])
        self.assertTrue(passed)
        self.assertGreater(min(pivots), 0)

    def test_prefix_master_permits_gain_above_one(self):
        c0, cn, qb = self.maps()
        args = (c0, cn, matrix_identity(18), matrix_identity(18),
                Interval.point(2), [(Interval.point(4), qb)])
        self.assertEqual(len(M.motion_master(*args, prefix=True)), 21)
        with self.assertRaises(ValueError):
            M.motion_master(*args)
        with self.assertRaises(ValueError):
            M.motion_master(c0[:17], cn[:17], matrix_identity(17), matrix_identity(17),
                            Interval.point(.5), [(Interval.point(4), qb)])

    def test_contract_requires_actual_active_filter_and_boundary_coverage(self):
        d = contract()
        self.assertEqual(M.validate(d), [])
        self.assertEqual(d["executed_active_filter_dimension"], 21)
        self.assertEqual(d["motion_storage_dimension"], 18)
        self.assertFalse(d["held_H18_certificate_transferred_to_active_filter"])
        self.assertFalse(d["existing_P3_automatically_covers_projection_boundary"])
        self.assertFalse(d["full_21_state_nonlinear_contraction_required"])
        self.assertFalse(d["sensor_noise_zero_implies_model_mismatch_zero"])
        self.assertFalse(d["zero_error_floor_required"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P5_MOTION_MAY_START"])
        for coverage in d["response_mode_coverage"].values():
            self.assertEqual(set(coverage), set(M.UNION.BRANCHES))
            self.assertTrue(all(c == {"H18": False, "A21": False} for c in coverage.values()))

    def test_symbolic_lemmas_cannot_forge_gains_or_deployment(self):
        d = contract()
        for field, value in (("P4_MOTION_PASS", True), ("P5_MOTION_MAY_START", True),
                             ("conditional_true_bias_bound_mps2", .4),
                             ("certified_motion_gains", {"rho": .9}),
                             ("sensor_noise_zero_implies_model_mismatch_zero", True),
                             ("empirical_accuracy_is_a_certified_theorem_constant", True)):
            bad = deepcopy(d)
            bad[field] = value
            self.assertTrue(M.validate(bad), field)


if __name__ == "__main__":
    unittest.main()
