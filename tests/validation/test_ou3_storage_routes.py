from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_storage_routes as R  # noqa: E402


class StorageRouteTests(unittest.TestCase):
    def test_congruent_coordinate_change_preserves_complete_word_ratio(self):
        rng = np.random.default_rng(916)
        a = rng.normal(size=(18, 18))/5
        p = rng.normal(size=(18, 18))
        q = rng.normal(size=(18, 18))
        m0, mn = p@p.T+np.eye(18), q@q.T+np.eye(18)
        units = np.diag(R.physical_scales(9.80665, 3.)[:18])
        expected, direction = R.largest_ratio(a, m0, mn)
        got, _ = R.largest_ratio(np.linalg.solve(units, a@units),
                                  units@m0@units, units@mn@units)
        self.assertAlmostEqual(expected, got, places=10)
        self.assertAlmostEqual(float(direction@m0@direction), 1., places=12)

    def example(self):
        rng = np.random.default_rng(917)
        steps, transitions, responses, measurements = [], [], [], []
        t, r = np.eye(21), np.zeros(21)
        for i in range(4):
            a = .75*np.eye(21)+rng.normal(scale=.02, size=(21, 21))
            b, u = rng.normal(scale=.03, size=(21, 3)), rng.normal(size=3)
            steps.append({"A": a, "B": b, "u": u, "index": i,
                          "stage": "projection", "kind": "accelerometer"})
            t, r = a@t, a@r+b@u
            transitions.append(t.copy())
            responses.append(r.copy())
            measurements.append({"H": rng.normal(size=(3, 21)),
                                 "R": np.diag([.1, .3, .7]), "kind": "accelerometer"})
        return steps, transitions, responses, measurements

    def test_signed_information_remainder_and_separate_bias_budget(self):
        steps, ts, rs, measurements = self.example()
        m0, metrics = np.eye(18), [np.eye(18)]*len(steps)
        centered = R.evaluate(steps, ts, rs, m0, metrics, measurements)
        result = R.transported_information(steps, ts, metrics, measurements, m0, centered)
        loss = np.eye(18)-ts[-1][:18, :18].T@ts[-1][:18, :18]
        np.testing.assert_allclose(result["signed_complete_word_loss_eigenvalues"],
                                   np.linalg.eigvalsh(loss), atol=2e-15)
        self.assertLess(result["matrix_identity_max_defect"], 1e-12)
        self.assertFalse(result["P4_PASS"])
        w = result["endpoint_witness"]
        info = sum(float(x) for x in w["transported_information_80_digit"].values())
        self.assertAlmostEqual(info+float(w["signed_remainder_80_digit"]),
                               1-float(w["ratio_80_digit"]), places=12)
        x, bias, amplitude = np.arange(18)/20, np.array([.2, -.1, .4]), -.7
        actual = ts[-1][:18]@np.r_[x, bias]+rs[-1][:18]*amplitude
        bound = (np.sqrt(centered["endpoint_ratio"])*np.linalg.norm(x)
                 + np.sqrt(centered["initial_bias_gain_squared_endpoint"])*np.linalg.norm(bias)
                 + np.sqrt(centered["particular_response_energy_endpoint"])*abs(amplitude))
        self.assertLessEqual(np.linalg.norm(actual), bound)

    def test_prefix_information_never_uses_future_measurements(self):
        steps, _, _, measurements = self.example()
        direction = np.arange(18)/20
        one = R.witness(steps, np.eye(18), [np.eye(18)]*4, direction, 1, measurements)
        measurements[2]["H"] *= 1e9
        two = R.witness(steps, np.eye(18), [np.eye(18)]*4, direction, 1, measurements)
        self.assertEqual(one, two)


if __name__ == "__main__":
    unittest.main()
