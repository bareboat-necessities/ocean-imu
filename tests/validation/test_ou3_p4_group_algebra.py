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
        # Both are independent outward enclosures of the same real quantity.
        self.assertLessEqual(max(direct.lo, identity.lo), min(direct.hi, identity.hi))

    def test_proof_module_does_not_use_linearized_attitude_injection(self):
        text = (ROOT / "tools" / "ou3_p4_group_algebra.py").read_text(encoding="utf-8")
        self.assertIn("deployed_injection_rotation", text)
        self.assertIn("_series_quaternion", text)
        self.assertIn("_axis_angle_quaternion", text)
        self.assertNotIn("I + [dtheta]", text)
        self.assertNotIn("math.sin(", text)
        self.assertNotIn("math.cos(", text)


if __name__ == "__main__":
    unittest.main()
