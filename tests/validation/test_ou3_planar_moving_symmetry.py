import unittest
from tools.stability.ou3_theorem.planar_moving_symmetry import certificate, mean_parity_residuals

class PlanarMovingSymmetryTests(unittest.TestCase):
    def test_parity_identities(self):
        c=certificate()
        self.assertLess(c["sampled_identity_worst_abs"],1e-12)
        self.assertFalse(c["all_time_magnetic_service_verified"])
    def test_exact_mag_residual_on_canonical_branch(self):
        q=mean_parity_residuals(.017,.002,-.001,[0.,0.],[0.,0.])
        self.assertLess(abs(q["mag_exact_residual_norm"]),1e-12)

if __name__=="__main__": unittest.main()
