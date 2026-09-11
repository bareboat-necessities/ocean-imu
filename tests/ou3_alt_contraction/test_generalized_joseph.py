import unittest
from fractions import Fraction as F
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import generalized_joseph as G

class GeneralizedJosephTests(unittest.TestCase):
    def test_masked_exact_regression(self):
        d=G.build(); self.assertEqual(G.validate(d),[])
        self.assertTrue(d['masked_regression_L_differs_from_physical_H'])
        self.assertTrue(d['masked_regression_C_differs_from_physical_R'])
        self.assertTrue(d['exact_masked_Joseph_signed_identity_closed'])
        self.assertTrue(d['exact_masked_Joseph_reset_signed_identity_closed'])
        self.assertFalse(d['source_uniform_complete_word_dissipation_closed_here'])

    def test_unmasked_limit_recovers_standard_information_form(self):
        d=G.build()
        self.assertTrue(d['unmasked_limit_recovers_L_equals_H'])
        self.assertTrue(d['unmasked_limit_recovers_C_equals_R'])
        self.assertFalse(d['unqualified_N_equals_P_Htranspose_required'])

    def test_descriptor_form_matches_inverse_form(self):
        m=G.masked_regression();s=m['signed'];d=m['descriptor']
        self.assertEqual(d['joseph_signed'],s['joseph_signed'])
        self.assertEqual(d['reset_reduced_signed'],s['reset_reduced_signed'])
        self.assertTrue(all(x==0 for x in d['S_u_minus_q']))
        self.assertTrue(all(x==0 for x in d['C_v_minus_xi']))
        self.assertTrue(all(x==0 for x in d['C_w_minus_Lb']))

    def test_arbitrary_actual_numerator_identity(self):
        P=[[F(4),F(1,3)],[F(1,3),F(3)]]
        S=[[F(5),F(1,4),F(0)],[F(1,4),F(6),F(1,5)],[F(0),F(1,5),F(7)]]
        N=[[F(1,3),F(-1,7),F(1,11)],[F(0),F(1,5),F(-1,13)]]
        t=G.signed_terms(P,S,N,[F(2,5),F(-1,6)],[F(1,4),F(-2,9),F(1,8)],[F(1,40),F(-1,50)])
        self.assertEqual(t['joseph_identity_residual'],0)
        self.assertEqual(t['reset_identity_residual'],0)

    def test_full_state_still_has_three_dimensional_core(self):
        P=G.eye(5);S=[[F(4),F(0),F(0)],[F(0),F(5),F(0)],[F(0),F(0),F(6)]]
        N=[[F(1,10),F(0),F(0)],[F(0),F(1,10),F(0)],[F(0),F(0),F(1,10)],[F(1,20),F(1,30),F(0)],[F(0),F(1,40),F(1,50)]]
        o=G.objects(P,S,N)
        self.assertEqual(G.shape(o['C']),(3,3));self.assertEqual(G.shape(o['L']),(3,5));self.assertEqual(G.shape(o['Jp']),(5,5))

if __name__=='__main__':unittest.main()
