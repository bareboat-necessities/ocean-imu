from pathlib import Path
import math
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_accelerometer_corotated_aw as C


def dot(a, b):
    return sum(float(x) * float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mv(A, x):
    return [sum(A[i][j] * x[j] for j in range(3)) for i in range(3)]


def rotation(axis, theta):
    n = norm(axis)
    u = [x / n for x in axis]
    x, y, z = u
    c = math.cos(theta)
    s = math.sin(theta)
    Cc = 1.0 - c
    return [
        [c + x*x*Cc, x*y*Cc - z*s, x*z*Cc + y*s],
        [y*x*Cc + z*s, c + y*y*Cc, y*z*Cc - x*s],
        [z*x*Cc - y*s, z*y*Cc + x*s, c + z*z*Cc],
    ]


class Ou3P4AccelerometerCorotatedAwTests(unittest.TestCase):
    def test_exact_residual_decomposition_and_aw_isometry(self):
        Rhat = rotation([0.3, -0.4, 0.5], 0.37)
        axis = [0.2, 0.3, 0.4]
        theta = 0.70
        E = rotation(axis, theta)
        Rtrue = mm(E, Rhat)

        ahat = [0.8, -1.1, 0.4]
        da = [0.7, -0.2, 0.5]
        bhat = [0.03, -0.04, 0.01]
        db = [-0.02, 0.01, 0.04]
        g = [0.0, 0.0, 9.80665]

        fhat = mv(Rhat, sub(ahat, g))
        predicted = add(fhat, bhat)
        truth = add(mv(Rtrue, sub(add(ahat, da), g)), add(bhat, db))
        residual_direct = sub(truth, predicted)

        Q = mm(mm(transpose(Rhat), E), Rhat)
        uaw = mv(Q, da)
        rotation_residual = sub(mv(E, fhat), fhat)
        residual_corotated = add(add(rotation_residual, mv(Rhat, uaw)), db)

        for x, y in zip(residual_direct, residual_corotated):
            self.assertAlmostEqual(x, y, places=12)
        self.assertAlmostEqual(norm(uaw), norm(da), places=12)

        QtQ = mm(transpose(Q), Q)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(QtQ[i][j], 1.0 if i == j else 0.0, places=12)

    def test_only_pure_rotation_remains_in_eta(self):
        axis = [0.2, 0.3, 0.4]
        n = norm(axis)
        u = [x / n for x in axis]
        theta = 0.70
        E = rotation(axis, theta)
        fhat = [1.2, -2.3, 9.0]
        cayley = [2.0 * math.tan(theta / 2.0) * x for x in u]

        yrot = sub(mv(E, fhat), fhat)
        hrot = cross(cayley, fhat)
        eta = sub(yrot, hrot)

        # Exact Cayley identities used by the signed Joseph route.
        self.assertAlmostEqual(dot(yrot, eta), 0.0, places=12)
        self.assertAlmostEqual(
            dot(eta, eta) / dot(yrot, yrot),
            dot(cayley, cayley) / 4.0,
            places=12,
        )

        # Adding arbitrary co-rotated a_w and accelerometer-bias errors to both
        # exact and tangent residuals leaves eta unchanged bit-for-formula.
        linear = [0.6, -0.1, 0.3]
        exact_full = add(yrot, linear)
        tangent_full = add(hrot, linear)
        eta_full = sub(exact_full, tangent_full)
        for x, y in zip(eta_full, eta):
            self.assertAlmostEqual(x, y, places=12)

    def test_source_bound_primitive_does_not_promote_P4(self):
        d = C.build()
        self.assertEqual([], C.validate(d))
        self.assertTrue(d["aw_error_exactly_linear_in_accelerometer_operation_coordinate"])
        self.assertEqual(0.0, d["latent_aw_nonlinear_eta_coefficient"])
        self.assertEqual(0.0, d["accelerometer_bias_nonlinear_eta_coefficient"])
        self.assertFalse(d["complete_H18_A21_word_established_here"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertGreaterEqual(d["outer_angle_rad"], 0.80)


if __name__ == "__main__":
    unittest.main()
