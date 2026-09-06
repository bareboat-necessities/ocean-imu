from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TOOLS / "stability"))

import ou3_p4_magnetometer_radial_joseph as M


def dot(a, b):
    return sum(float(x) * float(y) for x, y in zip(a, b))


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, s):
    return [s * x for x in a]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def transpose(A):
    return [list(row) for row in zip(*A)]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mv(A, x):
    return [sum(A[i][j] * x[j] for j in range(3)) for i in range(3)]


def skew(v):
    x, y, z = v
    return [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]]


def inverse3(A):
    a,b,c = A[0]
    d,e,f = A[1]
    g,h,i = A[2]
    C00=e*i-f*h; C01=-(d*i-f*g); C02=d*h-e*g
    C10=-(b*i-c*h); C11=a*i-c*g; C12=-(a*h-b*g)
    C20=b*f-c*e; C21=-(a*f-c*d); C22=a*e-b*d
    det=a*C00+b*C01+c*C02
    if abs(det) < 1e-15:
        raise ValueError("singular")
    # inverse is cofactor transpose / det
    return [[C00/det,C10/det,C20/det],[C01/det,C11/det,C21/det],[C02/det,C12/det,C22/det]]


def quad(x, A):
    return dot(x, mv(A, x))


class Ou3P4MagnetometerRadialJosephTests(unittest.TestCase):
    def test_radial_energy_cancels_exactly_in_example(self):
        v = [2.0, -1.0, 3.0]
        H = [[-x for x in row] for row in skew(v)]  # -[v]x
        P = [[2.0, 0.3, -0.2], [0.3, 1.5, 0.1], [-0.2, 0.1, 1.2]]
        r = 0.7
        S = mm(mm(H, P), transpose(H))
        for k in range(3):
            S[k][k] += r
        Sinv = inverse3(S)

        Sv = mv(S, v)
        for x, y in zip(Sv, scale(v, r)):
            self.assertAlmostEqual(x, y, places=11)

        y = [0.4, -0.7, 0.9]
        alpha = dot(y, v) / dot(v, v)
        radial = scale(v, alpha)
        tangent = sub(y, radial)
        c = [0.12, -0.08, 0.05]
        hlin = mv(H, c)
        eta = sub(y, hlin)
        eta_t = sub(tangent, hlin)

        full = quad(y, Sinv) - dot(eta, eta) / r
        reduced = quad(tangent, Sinv) - dot(eta_t, eta_t) / r
        self.assertAlmostEqual(full, reduced, places=11)

    def test_effective_tangent_coordinate_reconstructs_tangent_residual(self):
        v = [2.0, -1.0, 3.0]
        H = [[-x for x in row] for row in skew(v)]
        y = [0.4, -0.7, 0.9]
        alpha = dot(y, v) / dot(v, v)
        tangent = sub(y, scale(v, alpha))
        d = scale(mv(transpose(H), y), 1.0 / dot(v, v))
        Hd = mv(H, d)
        for x, z in zip(Hd, tangent):
            self.assertAlmostEqual(x, z, places=12)

    def test_source_bound_primitive_does_not_promote_P4(self):
        d = M.build()
        self.assertEqual([], M.validate(d))
        self.assertTrue(d["radial_Joseph_energy_cancellation_exact"])
        self.assertFalse(d["standalone_radial_eta_penalty_used"])
        self.assertTrue(d["directional_form_retained_until_word_scalarization"])
        self.assertFalse(d["complete_H18_A21_word_established_here"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
