from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.matrix_certificates import add, is_psd, ldlt
from tools.stability.ou3_theorem.root_covariance_certificate import certificate, process_floors


class JointRootCovarianceTests(unittest.TestCase):
    def test_convex_combination_preserves_carried_cross_covariance(self):
        p=[[F(2),F(1)],[F(1),F(2)]]
        x=[[F(3,2),F(0)],[F(0),F(0)]]
        y=[[F(0),F(0)],[F(0),F(3,2)]]
        self.assertTrue(is_psd(add(p,x,F(-1))))
        self.assertTrue(is_psd(add(p,y,F(-1))))
        self.assertFalse(is_psd(add(p,add(x,y),F(-1))))
        self.assertTrue(is_psd(add(p,add(x,y),F(-1,2))))

    def test_joint_matrix_is_strictly_positive_in_all_21_coordinates(self):
        r=certificate(); c=[[F(v) for v in row] for row in r['LIN_covariance_lower_matrix']]
        def q(record):return F(int(record['numerator']),int(record['denominator']))
        ag=q(r['AG_covariance_diagonal_lower']);ba=q(r['BA_covariance_diagonal_lower'])
        full=[[F(0) for _ in range(21)] for _ in range(21)]
        for i in range(6):full[i][i]=ag
        for i in range(3):full[18+i][18+i]=ba
        for axis in range(3):
            for i in range(4):
                for j in range(4):full[6+3*i+axis][6+3*j+axis]=c[i][j]
        _,d=ldlt(full)
        self.assertTrue(all(v>0 for v in d))
        self.assertEqual(ag,F('4.9991e-10'))
        self.assertGreater(ba,F('4.99999e-10'))
        self.assertFalse(r['constructive_full_A21_mu_rho_enclosure'])
        self.assertFalse(r['uniform_covariance_upper_bound_verified'])
        self.assertFalse(r['float32_covariance_factor_verified'])

    def test_source_gyro_floor_uses_deployed_density_not_old_simulation_value(self):
        source=(ROOT/'sensors/full_marine_ins/atomS3R_ins_kalman_ou3/atomS3R_ins_kalman_ou3.ino').read_text()
        self.assertIn('gyr_sigma_ref_rps  = 0.00135f',source)
        ag,ba=process_floors()
        self.assertGreater(ag,0)
        self.assertGreater(ba,0)


if __name__=='__main__':
    unittest.main()
