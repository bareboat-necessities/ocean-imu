from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_connected_motion as C  # noqa: E402
import ou3_p4_motion_gain as G  # noqa: E402


class MotionGainTests(unittest.TestCase):
    def test_gyro_factor_is_exact_finite_group_product(self):
        for a, b in (([0, 0, 0], [0, 0, 0]),
                     ([.001, -.002, .003], [1e-9, -2e-9, 3e-9]),
                     ([.001, 0, 0], [-.002, 0, 0]),
                     ([.02, -.03, .01], [.001, .004, -.002]),
                     ([.0099, 0, 0], [.0002, .001, 0])):
            a, b = np.array(a), np.array(b)
            got = G.gyro_secant(a, b)@b
            expected = C.cayley(C.deployed_rotation(a+b)@C.deployed_rotation(a).T)
            np.testing.assert_allclose(got, expected, atol=2e-16, rtol=1e-10)
        with self.assertRaises(ValueError):
            G.gyro_secant([.02, 0, 0], [0, 0, 0])

    def test_measurement_factor_includes_RS_forcing_finite_reset_and_projection(self):
        rng = np.random.default_rng(905)
        e = rng.normal(0, .1, 21)
        for kind in ("accelerometer", "magnetometer", "S_zero"):
            row = {"kind": kind, "R_hat": C.rotation([.1, -.03, .04]).ravel().tolist(),
                   "R_true": np.eye(3).ravel().tolist(), "x_hat": [0.]*21,
                   "linear_true": rng.normal(size=12).tolist(), "true_bias": [0.]*3,
                   "measured": [1., 2., 3.], "temperature_mean": [0.]*3,
                   "mag_reference": [20., 0., 40.], "gravity": 9.80665,
                   "projection_radius": .4, "R": np.diag([.1, .3, .7]).ravel().tolist(),
                   "K": rng.normal(0, .02, (21, 3)).ravel().tolist()}
            a, b, u, _, _ = G.measurement_factor(row, e)
            y, _, _ = C.residual_graph(row, e)
            expected = C.correct(e, C.mat(row, "K", 21, 3), y, [0.]*3, .4)
            np.testing.assert_allclose(a@e+b@u, expected, atol=2e-14)
            self.assertGreater(np.linalg.norm(b[:18]), 0)
            self.assertGreater(np.linalg.norm(a[18:, :18]), 0)
            if kind == "S_zero":
                np.testing.assert_array_equal(u, -np.array(row["linear_true"][6:9]))
            row["true_bias"] = [.01, 0, 0]
            with self.assertRaises(ValueError):
                G.measurement_factor(row, e)

    def test_bias_supply_elimination_keeps_dense_cross_terms(self):
        rng = np.random.default_rng(906)
        x = rng.normal(size=(21, 21))
        p = x@x.T+np.eye(21)
        expected = np.linalg.inv(np.linalg.inv(p)+np.diag([0.]*18+[3.]*3))
        np.testing.assert_allclose(G.supply_update(p, 3.), expected, atol=2e-13)

    def test_prediction_factor_keeps_physical_forcing_and_finite_gyro_bias(self):
        rng = np.random.default_rng(908)
        e = rng.normal(0, .05, 21)
        before = {"x_hat": (-e).tolist(), "R_true": C.rotation(e[:3]).ravel().tolist(),
                  "R_hat": np.eye(3).ravel().tolist(), "P": np.eye(21).ravel().tolist(),
                  "linear_true": [0.]*12, "true_bias": [0.]*3}
        dt, omega = .005, np.array([.3, -.1, .2])
        f = .99*np.eye(12)+rng.normal(0, .001, (12, 12))
        u, ug = rng.normal(0, .001, 12), np.array([.0001, -.0002, .0003])
        rs, rn = C.deployed_rotation(-dt*(omega-e[3:6])), C.deployed_rotation(-dt*omega)
        for active in (False, True):
            row = {"active": active, "tau_b": 5000., "F_LL": f.ravel().tolist(),
                   "omega_hat": omega.tolist(), "gyro_measured": (omega-e[3:6]).tolist(),
                   "linear_true": u.tolist(),
                   "R_true": (C.rotation(ug)@rs@C.rotation(e[:3])).ravel().tolist()}
            a, b, source, _, _ = G.prediction_factor(before, row, dt)
            expected = e.copy()
            expected[:3] = C.cayley(C.rotation(ug)@rs@C.rotation(e[:3])@rn.T)
            expected[6:18] = f@e[6:18]+u
            expected[18:] *= np.exp(-dt/5000.) if active else 1.
            np.testing.assert_allclose(a@e+b@source, expected, atol=2e-15)

    def example(self):
        rng = np.random.default_rng(907)
        steps = []
        for i in range(6):
            entrance = i in (0, 3)
            a = np.eye(21) if entrance else .8*np.eye(21)+rng.normal(0, .02, (21, 21))
            b = np.zeros((21, 0)) if entrance else rng.normal(0, .03, (21, 3))
            steps.append({"A": a, "B": b, "Uinv": np.eye(b.shape[1]),
                          "metric": np.diag(np.linspace(1, 2, 18)),
                          "bias_cost": entrance, "stage": "prediction_enter" if entrance else "projection",
                          "index": i, "kind": "S_zero"})
        return steps, np.eye(18)

    def test_recursive_gain_equals_full_augmented_matrix_with_corrected_bias(self):
        steps, m0 = self.example()
        dt, factor, gain = .005, .9, 20.
        size = 21+sum(s["B"].shape[1] for s in steps)
        state = np.zeros((21, size))
        state[:, :21] = np.eye(21)
        cost = np.zeros((size, size))
        cost[:18, :18] = factor*m0
        offset = 21
        expected = []
        for s in steps:
            if s["bias_cost"]:
                cost += gain*dt*state[18:].T@state[18:]
            else:
                n = s["B"].shape[1]
                state = s["A"]@state
                state[:, offset:offset+n] += s["B"]
                cost[offset:offset+n, offset:offset+n] = gain*np.linalg.inv(s["Uinv"])
                offset += n
            used = cost[:offset, :offset]
            chol = np.linalg.cholesky(s["metric"])
            reach = state[:18, :offset]@np.linalg.solve(used, state[:18, :offset].T)
            expected.append(np.linalg.eigvalsh(G.sym(chol.T@reach@chol))[-1])
        actual = G.gain_test(steps, m0, dt, factor, gain, retain=True)
        np.testing.assert_allclose(actual["ratios"], expected, atol=3e-15)
        witness = G.high_precision_direction(steps, m0, dt, factor, gain, actual)
        self.assertAlmostEqual(float(witness["ratio_80_digit"]), actual["endpoint"], places=12)
        self.assertAlmostEqual(float(witness["supply_cost_80_digit"]), 1., places=12)
        self.assertFalse(witness["BRMM_nonlinear_admissibility_of_maximizer_established"])

    def test_bad_supply_or_unattached_word_fails_closed(self):
        steps, m0 = self.example()
        with self.assertRaises(ValueError):
            G.gain_test(steps, m0, .005, 1., 0.)
        with self.assertRaises(ValueError):
            G.build_word({"dt": .005}, [], [], "A21")


if __name__ == "__main__":
    unittest.main()
