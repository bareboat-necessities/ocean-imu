from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import solve_triangular

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.factor_word import FactorWord, run_word
from tools.stability.ou3_theorem.factor_certificates import eliminate_nuisance, singular_floor
from tools.stability.ou3_theorem.matrix_certificates import identity
from tools.stability.ou3_theorem.word_energy import word_identity


class FactorWordTests(unittest.TestCase):
    def compare(self,p,events):
        exact=word_identity(p,events)
        with patch('numpy.linalg.inv',side_effect=AssertionError('dense inverse')):
            word=run_word(p,events)
        c=word.root_factor
        loss=c.T@np.array(exact['loss'],float)@c
        np.testing.assert_allclose(word.loss.T@word.loss,loss,rtol=2e-12,atol=2e-14)
        np.testing.assert_allclose(word.c@word.c.T,np.array(exact['end_covariance'],float),rtol=2e-12,atol=2e-14)
        np.testing.assert_allclose(word.t,np.array(exact['end_transport'],float)@c,rtol=2e-12,atol=2e-14)
        self.assertLess(word.snapshot()['identity_residual_fro'],1e-12)
        self.assertFalse(word.snapshot()['source_uniform_verified'])
        return word

    def test_interleaved_correlated_noise_singular_increment_and_nonorthogonal_reset(self):
        self.compare([[2,1],[1,3]],[
            {'kind':'prediction','F':[[1,F(1,5)],[0,1]],'Q':[[1,F(1,3)],[F(1,3),1]]},
            {'kind':'correction','H':[[1,2]],'R':[[F(3,2)]]},
            {'kind':'reset','G':[[1,F(1,10)],[-F(1,10),1]]},
            {'kind':'prediction','F':[[1,0],[0,F(9,10)]],'Q':[[0,0],[0,2]]},
            {'kind':'correction','H':[[1,0],[0,1]],'R':[[2,1],[1,2]]}])

    def test_full_21_state_cross_coupling_and_three_row_measurement(self):
        p=identity(21);p[0][20]=p[20][0]=F(1,3)
        f=identity(21);f[0][3]=F(1,200)
        h=[[F(0) for _ in range(21)] for _ in range(3)]
        h[0][0]=h[0][20]=h[1][1]=h[2][2]=F(1)
        q=[[F(v,100) for v in row] for row in identity(21)]
        w=self.compare(p,[{'kind':'prediction','F':f,'Q':q},
                          {'kind':'correction','H':h,'R':identity(3)}])
        self.assertEqual(w.snapshot()['max_measurement_rows'],3)
        self.assertEqual(w.loss.shape,(21,21))
        self.assertGreater(w.snapshot()['delta_diagnostic'],0)

    def test_singular_transition_with_sufficient_process_excitation(self):
        self.compare(identity(2),[{'kind':'prediction','F':[[1,0],[0,0]],'Q':identity(2)}])

    def test_cancellation_survives_factorization_and_rejected_event(self):
        w=run_word(identity(3),[
            {'kind':'correction','H':[[1,0,1],[0,1,0]],'R':identity(2)},
            {'kind':'correction','H':identity(3),'R':identity(3),'accepted':False}])
        self.assertLess(np.linalg.norm(w.loss@np.array([1,0,-1])),1e-14)
        self.assertEqual(w.snapshot()['delta_diagnostic'],0)
        self.assertEqual(w.counts['rejected'],1)

    def test_innovation_safety_boost_is_effective_noise(self):
        w=FactorWord(identity(2));w.correction([[1,1]],[[1]],innovation_increment=[[F(1,10)]])
        ref=self.compare(identity(2),[{'kind':'correction','H':[[1,1]],'R':[[F(11,10)]]}])
        np.testing.assert_allclose(w.loss.T@w.loss,ref.loss.T@ref.loss,atol=1e-14)
        with self.assertRaises(ValueError):
            w.correction([[1,1]],[[1]],innovation_increment=[[-1]])
        with self.assertRaises(ValueError):
            w.correction([[1,1]],[[1]],innovation_increment=identity(2))
        with self.assertRaises(ValueError):
            FactorWord([[1,1],[0,1]])  # a lower-only Cholesky would ignore the defect


class FactorCertificateTests(unittest.TestCase):
    def test_exact_singular_nuisance_cancellation(self):
        r=eliminate_nuisance([[1,0,1,1],[0,1,0,0]],[0,1])
        self.assertEqual(r['reduced_information'],[[0,0],[0,1]])
        self.assertEqual(r['minimizer_map'],[[-1,0],[0,0]])
        self.assertTrue(r['nuisance_gram_singular'])

    def test_positive_heading_schur_does_not_remove_nuisance_nullspace(self):
        r=eliminate_nuisance([[1,0,1,1],[0,1,0,0],[0,0,1,1]],[0,1])
        self.assertEqual(r['reduced_information'],[[F(1,2),0],[0,1]])
        self.assertEqual(r['full_factor_rank'],3)  # four state columns

    def test_rational_rank_has_no_small_singular_value_cutoff(self):
        eps=F(1,10**40)
        r=eliminate_nuisance([[1,1,1],[0,0,eps],[0,0,0]],[0])
        self.assertEqual(r['nuisance_rank'],2)
        self.assertEqual(r['reduced_information'],[[0]])

    def test_residual_floor_requires_full_factor_error_budget(self):
        b=np.array([[2.,1.],[0.,3.]])
        x=solve_triangular(b,np.eye(2))
        good=singular_floor(b,x,factor_error_norm_upper=F(1,100))
        self.assertTrue(good['positive_floor'])
        self.assertLess(float(F(good['sigma_lower'])),np.linalg.svd(b,compute_uv=False)[-1])
        self.assertFalse(good['whole_word_error_verified'])
        bad=singular_floor(b,x,factor_error_norm_upper=3)
        self.assertFalse(bad['positive_floor'])
        singular=singular_floor([[1,1],[0,0]],identity(2),factor_error_norm_upper=0)
        self.assertFalse(singular['positive_floor'])


if __name__=='__main__':
    unittest.main()
