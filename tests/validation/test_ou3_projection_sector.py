"""Regression tests for the exact-real joint radial projection sector."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[2] / "tools/stability/ou3_projection_sector.py"
SPEC = importlib.util.spec_from_file_location("projection_sector", PATH)
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


class ProjectionSectorTests(unittest.TestCase):
    def test_inside_and_outside_jacobians_have_joint_gain_at_most_one(self):
        r = .4
        for x in ([.1, .2, 0.], [.8, 0., 0.], [.5, -.7, .3]):
            j = P.projection_jacobian(x, r)
            np.testing.assert_allclose(j, j.T, rtol=0, atol=2e-15)
            ev = np.linalg.eigvalsh(j)
            self.assertGreaterEqual(ev[0], -2e-15)
            self.assertLessEqual(ev[-1], 1+2e-15)
            self.assertLessEqual(P.stacked_sector_ratio(j), 1+2e-15)

    def test_boundary_clarke_family_has_same_sector(self):
        tangent = np.diag([0., 1., 1.])
        for theta in np.linspace(0, 1, 21):
            j = theta*np.eye(3)+(1-theta)*tangent
            self.assertLessEqual(P.stacked_sector_ratio(j), 1+2e-15)

    def test_joint_finite_difference_sector_crosses_saturation_boundary(self):
        # These pairs deliberately place beta-e on different sides/directions
        # of the projection boundary.  This is a regression, while the module
        # docstring contains the all-points generalized-Jacobian proof.
        r = .4
        pairs = [
            (np.array([.05, 0, 0]), np.zeros(3),
             np.array([-.7, .2, 0]), np.array([.15, -.1, 0])),
            (np.array([.8, -.2, .1]), np.array([.1, .2, -.1]),
             np.array([-.5, .6, -.4]), np.array([-.2, .1, .3])),
            (np.array([.4, 0, 0]), np.zeros(3),
             np.array([.400001, 0, 0]), np.zeros(3)),
        ]
        for e1,b1,e2,b2 in pairs:
            f1=P.bias_error_projection(e1,b1,r)
            f2=P.bias_error_projection(e2,b2,r)
            lhs=float(np.dot(f1-f2,f1-f2))
            rhs=float(np.dot(e1-e2,e1-e2)+np.dot(b1-b2,b1-b2))
            self.assertLessEqual(lhs, rhs*(1+2e-14)+1e-15)

    def test_fixed_physical_bias_is_nonexpansive_in_error(self):
        r=.4
        beta=np.array([.3,-.2,.1])
        rng=np.random.default_rng(7)
        for _ in range(200):
            e1,e2=rng.normal(size=(2,3))
            f1=P.bias_error_projection(e1,beta,r)
            f2=P.bias_error_projection(e2,beta,r)
            self.assertLessEqual(np.linalg.norm(f1-f2),
                                 np.linalg.norm(e1-e2)*(1+2e-14)+1e-15)

    def test_estimate_ball_invariant_in_all_regimes(self):
        r=.4
        for e,b in (([0,0,0],[0,0,0]), ([2,-3,4],[.2,-.1,.05]),
                    ([-.1,.2,.3],[1.,-1.,.5])):
            projected_error=P.bias_error_projection(e,b,r)
            estimate=np.asarray(b)-projected_error
            self.assertLessEqual(np.linalg.norm(estimate),r*(1+2e-15))

    def test_report_never_promotes_p4(self):
        d=P.build_report(.4)
        self.assertTrue(d["exact_real_operator_sector_closed"])
        self.assertTrue(d["saturated_branch_included_analytically"])
        self.assertTrue(d["boundary_clarke_branch_included_analytically"])
        self.assertFalse(d["fixed_multiplier_required"])
        for key in ("P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
            self.assertFalse(d[key])

    def test_input_validation(self):
        for r in (-1,np.nan,np.inf):
            with self.subTest(r=r), self.assertRaises(ValueError):
                P.project_ball([1,2,3],r)
        with self.assertRaises(ValueError):
            P.projection_jacobian([.4,0,0],.4)


if __name__ == "__main__":
    unittest.main()
