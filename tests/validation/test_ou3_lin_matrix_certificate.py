from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.lin_matrix_certificate import (
    SCALES, action_matrix, certificate, gram, gram_ceiling, translation_normalization,
)
from tools.stability.ou3_theorem.matrix_certificates import add, congruence, identity, is_psd, ldlt, matmul, transpose


class MatrixPathTests(unittest.TestCase):
    def test_exact_ldlt_reconstructs_positive_precision(self):
        a=action_matrix(); l,d=ldlt(a)
        diag=[[d[i] if i==j else F(0) for j in range(4)] for i in range(4)]
        self.assertEqual(congruence(diag,transpose(l)),a)
        r=certificate(); inv=[[F(v) for v in row] for row in r['covariance_lower_matrix']]
        self.assertEqual(matmul(a,inv),identity(4))

    def test_time_envelope_contains_endpoint_and_interior_gram_matrices(self):
        # Consistency check of the analytic monotone-entry/row-sum envelope.
        scales=tuple(map(F,SCALES))
        for order in range(5):
            upper=gram_ceiling(order,scales)
            for time in (F(16),F('16.003'),F('16.006')):
                self.assertTrue(is_psd(add(upper,gram(order,time,scales),F(-1))))

    def test_source_uniform_event_delay_normalization_is_exact_and_positive(self):
        r=certificate(); mu=r['normalization_raw_translation_lower']
        self.assertGreater(F(int(mu['numerator']),int(mu['denominator'])),F('3.59e-11'))
        self.assertFalse(r['full_state_coercivity_verified'])
        self.assertFalse(r['constructive_full_A21_mu_rho_enclosure'])

    def test_congruent_coordinate_change_does_not_manufacture_information(self):
        a=action_matrix(); scale=[F(2),F(3),F(5),F(7)]
        new_a=[[a[i][j]*scale[i]*scale[j] for j in range(4)] for i in range(4)]
        new_scales=[F(s)*d for s,d in zip(SCALES,scale)]
        old=translation_normalization(a); new=translation_normalization(new_a,new_scales)
        self.assertLess(abs(new/old-1),F('1e-20'))

    def test_ldlt_pivots_are_not_eigenvalue_floors(self):
        a=[[F(1),F(1)],[F(1),F(2)]]
        _,d=ldlt(a)
        self.assertEqual(d,[1,1])
        self.assertFalse(is_psd(add(a,identity(2),F(-1))))


if __name__=='__main__':
    unittest.main()
