"""Algebraic regressions; none of these tests promotes nonlinear P4."""
from __future__ import annotations

from decimal import Decimal
import importlib.util
from pathlib import Path
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[1] / "kalman_ou_iii/ou3_p4_bounded_bias_cocycle.py"
SPEC = importlib.util.spec_from_file_location("bounded_bias_cocycle", PATH)
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


def step(a=None, *, index=0, entrance=False, force=None, kind="none"):
    a = np.eye(21) if a is None else np.asarray(a, dtype=float)
    force = np.zeros(21) if force is None else np.asarray(force, dtype=float)
    return {"A": a, "B": force[:, None], "u": np.ones(1),
            "bias_cost": entrance, "index": index,
            "stage": "prediction_enter" if entrance else "projection",
            "kind": kind}


def stable_sample(index=0):
    a = np.eye(21)
    a[:18, :18] *= .5
    a[0, 18] = .2
    force = np.zeros(21)
    force[1] = .1
    return [step(index=index, entrance=True), step(a, index=index, force=force)]


class SampleLiftTests(unittest.TestCase):
    def test_full_sample_composition_preserves_internal_bias_feedback(self):
        to_bias, to_motion = np.eye(21), np.eye(21)
        to_bias[18, 0] = 2
        to_motion[0, 18] = 3
        samples = C.sample_lifts([step(entrance=True), step(to_bias), step(to_motion)])
        actual = C.motion_scan(samples, .01)[-1]["transition"]
        self.assertEqual(actual[0, 0], 7)
        wrong = to_motion[:18, :18] @ to_bias[:18, :18]
        self.assertEqual(wrong[0, 0], 1)

    def test_between_sample_selection_is_not_full_word_principal_block(self):
        to_bias, to_motion = np.eye(21), np.eye(21)
        to_bias[18, 0] = 2
        to_motion[0, 18] = 3
        samples = C.sample_lifts([step(entrance=True), step(to_bias),
                                  step(index=1, entrance=True), step(to_motion, index=1)])
        scan = C.motion_scan(samples, .01)
        self.assertEqual(scan[-1]["transition"][0, 0], 1)
        self.assertEqual((to_motion @ to_bias)[0, 0], 7)
        initial = np.zeros(21)
        initial[0] = 1
        self.assertLess(C.reconstruction_defect(samples, initial), 1e-14)

    def test_held_cascade_agrees_with_full_word(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        full = samples[1]["full"] @ samples[0]["full"]
        np.testing.assert_allclose(C.motion_scan(samples, .01)[-1]["transition"],
                                   full[:18, :18])

    def test_bias_cost_once_per_entrance_matches_explicit_columns(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        scan = C.motion_scan(samples, .25)
        a, g = samples[0]["full"][:18, :18], samples[0]["full"][:18, 18:]
        columns = np.concatenate((a @ g, g), axis=1) / np.sqrt(.25)
        np.testing.assert_allclose(scan[-1]["gramian"], columns @ columns.T)
        np.testing.assert_allclose(scan[2]["gramian"], scan[1]["gramian"])

    def test_same_template_is_one_correlated_column(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        response = C.motion_scan(samples, .01)[-1]["response"]
        expected = samples[1]["full"][:18, :18] @ samples[0]["forcing"][:18]
        expected += samples[1]["forcing"][:18]
        np.testing.assert_allclose(response, expected)
        self.assertAlmostEqual(response[1]**2, .0225)
        self.assertNotAlmostEqual(response[1]**2, .1**2 + .05**2)

    def test_every_completed_prefix_is_retained(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        self.assertEqual(len(C.motion_scan(samples, .01)), 4)

    def test_word_and_entrance_validation(self):
        for steps in ([], [step()], [step(entrance=True), step(index=2, entrance=True)],
                      [step(np.eye(18), entrance=True)],
                      [step(2*np.eye(21), entrance=True)]):
            with self.subTest(steps=len(steps)), self.assertRaises(ValueError):
                C.sample_lifts(steps)
        invalid = np.eye(21)
        invalid[0, 0] = np.nan
        with self.assertRaises(ValueError):
            C.sample_lifts([step(entrance=True), step(invalid)])
        for dt in (0, -1, np.nan, np.inf):
            with self.subTest(dt=dt), self.assertRaises(ValueError):
                C.motion_scan(C.sample_lifts(stable_sample()), dt)


class StorageTests(unittest.TestCase):
    def test_constructive_metric_supports_defective_stable_matrix(self):
        a = np.array([[.5, 2.], [0., .5]])
        metric, report = C.stein_metric(a, np.ones(2))
        q = report["stein_q"]
        np.testing.assert_allclose(metric-(a/q).T @ metric @ (a/q), np.eye(2), atol=1e-12)
        self.assertLess(report["complete_word_ratio"], 1)
        self.assertGreater(np.linalg.eigvalsh(metric)[0], 0)

    def test_unstable_and_nonfinite_inputs_fail(self):
        for a in (np.eye(2), 1.01*np.eye(2), np.full((2, 2), np.nan)):
            with self.subTest(a=a), self.assertRaises(ValueError):
                C.stein_metric(a, np.ones(2))
        with self.assertRaises(ValueError):
            C.stein_metric(.5*np.eye(2), np.zeros(2))

    def test_per_word_schur_does_not_imply_consecutive_stability(self):
        u = np.array([[.5, 2.], [0., .5]])
        v = u.T
        self.assertLess(max(abs(np.linalg.eigvals(u))), 1)
        self.assertLess(max(abs(np.linalg.eigvals(v))), 1)
        self.assertGreater(max(abs(np.linalg.eigvals(v @ u))), 1)

    def test_separate_supply_matches_dense_quadratic_master(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        end = C.motion_scan(samples, .25)[-1]
        metric, report = C.stein_metric(end["transition"], np.ones(18))
        supply = C.supply_constants(end, metric, report["complete_word_ratio"])
        a, g = samples[0]["full"][:18, :18], samples[0]["full"][:18, 18:]
        b = np.concatenate((a @ g, g), axis=1) / np.sqrt(.25)
        graph = np.column_stack((end["transition"], b, end["response"]))
        cost = np.zeros((25, 25))
        cost[:18, :18] = supply["rho"]*metric
        cost[18:24, 18:24] = supply["gamma_bias"]*np.eye(6)
        cost[24, 24] = supply["gamma_template"]
        self.assertLess(np.linalg.eigvalsh(graph.T @ metric @ graph-cost)[-1], 1e-12)
        self.assertLess(supply["normalized_supply_ratio"], 1)

    def test_zero_supply_channels_are_omitted(self):
        end = {"transition": .5*np.eye(18), "gramian": np.zeros((18, 18)),
               "response": np.zeros(18)}
        supply = C.supply_constants(end, np.eye(18), .25)
        self.assertEqual(supply["gamma_bias"], 0)
        self.assertEqual(supply["gamma_template"], 0)
        self.assertEqual(supply["normalized_supply_ratio"], 0)

    def test_decimal_direction_and_signed_costs(self):
        samples = C.sample_lifts(stable_sample() + stable_sample(1))
        x = np.zeros(18)
        x[0] = 1
        answer = C.decimal_direction(samples, x, np.eye(18))
        self.assertEqual(Decimal(answer["ratio"]), Decimal("0.0625"))
        total = sum(Decimal(v) for v in answer["signed_operation_costs"].values())
        self.assertEqual(total, Decimal(answer["ratio"])-1)
        prefix = C.decimal_direction(samples, x, np.eye(18), stop=1)
        self.assertEqual(Decimal(prefix["ratio"]), Decimal("0.25"))

    def test_audit_never_promotes_p4_or_p5(self):
        report = C.audit_steps(stable_sample() + stable_sample(1), np.zeros(21), .25,
                               np.ones(18))
        self.assertTrue(report["point_complete_word_feasible"])
        for name in ("P4_PASS", "P4_MOTION_PASS", "P5_MAY_START"):
            self.assertFalse(report[name])
        self.assertIn("endpoint_direction", report)
        self.assertIn("prefix_direction", report)


class ProjectionTests(unittest.TestCase):
    def test_projection_maps_outside_points_to_ball(self):
        for value in (np.zeros(3), [.4, 0, 0], [.400126, 0, 0], [2., -3., 4.]):
            out = C.project_ball(value, .4)
            self.assertLessEqual(np.linalg.norm(out), .4*(1+2e-15))

    def test_true_bias_may_be_outside_estimate_ball(self):
        true = np.array([1., -2., 3.])
        estimate = C.project_ball([-10., 20., -30.], .4)
        self.assertLessEqual(np.linalg.norm(true-estimate), np.linalg.norm(true)+.4+1e-14)

    def test_frozen_interior_factor_is_not_boundary_projection(self):
        matrix = np.eye(3)*1.0003146
        boundary = np.array([.4, 0, 0])
        self.assertGreater(np.linalg.norm(matrix @ boundary), .4)
        self.assertLessEqual(np.linalg.norm(C.project_ball(matrix @ boundary, .4)), .4)

    def test_projection_input_validation_and_overflow_avoidance(self):
        self.assertLessEqual(np.linalg.norm(C.project_ball([1e308, -1e308, 1e308], .4)), .4+1e-15)
        np.testing.assert_array_equal(C.project_ball([1., 2., 3.], 0), np.zeros(3))
        for value, radius in (([np.nan], .4), ([1.], -1), ([1.], np.inf)):
            with self.subTest(value=value), self.assertRaises(ValueError):
                C.project_ball(value, radius)


if __name__ == "__main__":
    unittest.main()
