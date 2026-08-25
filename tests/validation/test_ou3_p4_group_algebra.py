from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, matrix_mul
import ou3_p4_group_algebra as G


def p(x):
    return Interval.point(float(x))


class Ou3P4GroupAlgebraTests(unittest.TestCase):
    def test_deployed_series_branch_contains_source_formula_rotation(self):
        d = [0.005, -0.002, 0.001]
        R = G.deployed_injection_rotation(G.point_vector(d))
        t2 = sum(x*x for x in d)
        t4 = t2*t2
        w = 1.0 - t2/8.0 + t4/384.0
        k = 0.5 - t2/48.0 + t4/3840.0
        q = [w, k*d[0], k*d[1], k*d[2]]
        n = math.sqrt(sum(x*x for x in q))
        q = [x/n for x in q]
        w,x,y,z = q
        exact = [
            [1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
            [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
            [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)],
        ]
        for i in range(3):
            for j in range(3):
                self.assertTrue(R[i][j].contains(exact[i][j]))

    def test_axis_angle_branch_contains_rodrigues_rotation(self):
        d = [0.2, 0.0, 0.0]
        R = G.deployed_injection_rotation(G.point_vector(d))
        exact = [
            [1.0, 0.0, 0.0],
            [0.0, math.cos(0.2), -math.sin(0.2)],
            [0.0, math.sin(0.2), math.cos(0.2)],
        ]
        for i in range(3):
            for j in range(3):
                self.assertTrue(R[i][j].contains(exact[i][j]))

    def test_branch_straddling_box_hulls_both_source_branches(self):
        d = [Interval(0.0099, 0.0101), p(0.0), p(0.0)]
        R = G.deployed_injection_rotation(d)
        self.assertLessEqual(R[1][1].lo, math.cos(0.0101))
        self.assertGreaterEqual(R[1][1].hi, math.cos(0.0099))

    def test_group_energy_is_exact_trace_energy(self):
        R = G.deployed_injection_rotation(G.point_vector([0.3,0.0,0.0]))
        V = G.group_energy(R)
        self.assertTrue(V.contains(1.0-math.cos(0.3)))

    def test_exact_energy_identity_matches_direct_group_product(self):
        Re = G.rodrigues_rotation(G.point_vector([0.1,-0.03,0.02]))
        d = G.point_vector([-0.02,0.01,-0.005])
        before = G.group_energy(Re)
        after = G.group_energy(matrix_mul(G.rodrigues_rotation(d), Re))
        direct = after - before
        identity = G.exact_energy_change_identity(Re, d)
        self.assertLessEqual(max(direct.lo, identity.lo), min(direct.hi, identity.hi))

    def test_cayley_coordinate_is_two_tan_half_angle(self):
        theta = 0.7
        R = G.rodrigues_rotation(G.point_vector([theta,0.0,0.0]))
        c = G.cayley_coordinate(R)
        exact = 2.0*math.tan(theta/2.0)
        self.assertTrue(c[0].contains(exact))
        self.assertTrue(c[1].contains(0.0))
        self.assertTrue(c[2].contains(0.0))

    def test_inverse_cayley_round_trip_contains_rotation(self):
        c = G.point_vector([0.2,-0.05,0.03])
        R = G.rotation_from_cayley(c)
        cc = G.cayley_coordinate(R)
        for i, exact in enumerate((0.2,-0.05,0.03)):
            self.assertTrue(cc[i].contains(exact))

    def test_left_cayley_composition_matches_matrix_product(self):
        ca = G.point_vector([0.15,-0.02,0.01])
        cb = G.point_vector([-0.04,0.08,0.03])
        cc = G.cayley_compose_left(ca, cb)
        Rprod = matrix_mul(G.rotation_from_cayley(ca), G.rotation_from_cayley(cb))
        cref = G.cayley_coordinate(Rprod)
        for i in range(3):
            self.assertLessEqual(max(cc[i].lo, cref[i].lo), min(cc[i].hi, cref[i].hi))

    def test_deployed_injection_cayley_matches_deployed_rotation(self):
        for d in ([0.005,-0.002,0.001], [0.2,0.03,-0.01]):
            with self.subTest(d=d):
                D = G.point_vector(d)
                c = G.deployed_injection_cayley(D)
                R = G.deployed_injection_rotation(D)
                cref = G.cayley_coordinate(R)
                for i in range(3):
                    self.assertLessEqual(max(c[i].lo, cref[i].lo), min(c[i].hi, cref[i].hi))

    def test_proof_module_does_not_use_linearized_attitude_injection(self):
        text = (ROOT / "tools" / "ou3_p4_group_algebra.py").read_text(encoding="utf-8")
        self.assertIn("deployed_injection_cayley", text)
        self.assertIn("cayley_compose_left", text)
        self.assertIn("_series_quaternion", text)
        self.assertIn("_axis_angle_quaternion", text)
        self.assertNotIn("I + [dtheta]", text)
        self.assertNotIn("math.sin(", text)
        self.assertNotIn("math.cos(", text)


if __name__ == "__main__":
    unittest.main()
