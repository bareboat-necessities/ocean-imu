from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.lin_path_certificate import (
    certificate, derivative, hermite_basis, small_x_source_defect, sqrt_floor,
)


class LinPathCertificateTests(unittest.TestCase):
    def test_endpoint_paths_have_exact_jets(self):
        for p, selected in zip(hermite_basis(), (2,1,0,3)):
            for d in range(4):
                q=derivative(p,d)
                self.assertEqual(q[0],0)
                self.assertEqual(sum(q),int(d==selected))

    def test_source_polynomial_relative_error_has_strict_margin(self):
        error,j,_,_=small_x_source_defect()
        self.assertGreater(error,0)
        self.assertLess(error,F('0.00002'))
        self.assertGreater(j,0)

    def test_rational_factor_covers_full_mesh_and_tau_domain(self):
        r=certificate()
        self.assertTrue(r['verified'])
        self.assertGreater(F(r['factor_floor']),F('0.00000119'))
        p=r['covariance_floor']
        floor=F(int(p['numerator']),int(p['denominator']))
        self.assertLessEqual(F(r['factor_floor'])**2,floor)
        self.assertTrue(r['arbitrary_piecewise_constant_tau'])
        self.assertTrue(r['all_acc_and_integral_corrections_included'])
        self.assertFalse(r['float32_covariance_factor_verified'])
        self.assertFalse(r['constructive_full_A21_mu_rho_enclosure'])
        self.assertFalse(r['theorem_closed'])

    def test_square_root_export_is_outward(self):
        for q in (F(2),F(1,3),F(1,10**30)):
            x=sqrt_floor(q)
            self.assertLessEqual(x*x,q)
            self.assertGreater((x+F(1,10**24))**2,q)


if __name__=='__main__':
    unittest.main()
