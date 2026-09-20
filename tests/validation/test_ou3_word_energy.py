from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.matrix_certificates import congruence, identity, is_psd
from tools.stability.ou3_theorem.word_energy import (
    full_loss_margin, restricted_service_counterexample, word_identity,
)


class WordEnergyTests(unittest.TestCase):
    def test_restricted_service_cannot_be_lifted_to_full_heading_information(self):
        r=restricted_service_counterexample()
        self.assertTrue(r['verified'])
        self.assertEqual(r['restricted_heading_bias_loss'],[['1','0'],['0','1']])
        self.assertEqual(r['root_energy'],r['endpoint_energy'])
        self.assertFalse(r['full_state_strict_margin'])
        d=[[F(x) for x in row] for row in r['full_loss']]
        d[0][0]-=1
        self.assertFalse(is_psd(d))

    def test_complete_loss_telescopes_through_interleaved_operations(self):
        r=word_identity([[2,1],[1,3]],[
            {'kind':'prediction','F':[[1,F(1,5)],[0,1]],'Q':[[1,F(1,3)],[F(1,3),1]]},
            {'kind':'correction','H':[[1,2]],'R':[[F(3,2)]]},
            {'kind':'reset','G':[[1,F(1,10)],[-F(1,10),1]]},
            {'kind':'prediction','F':[[1,0],[0,F(9,10)]],'Q':[[0,0],[0,2]]},
            {'kind':'correction','H':[[1,0],[0,1]],'R':[[1,0],[0,2]]},
        ])
        self.assertTrue(r['algebra_verified'])
        self.assertTrue(full_loss_margin(r,F(1,100)))
        self.assertFalse(r['source_uniform_verified'])

    def test_null_direction_survives_even_with_positive_restricted_loss(self):
        r=word_identity(identity(2),[{'kind':'correction','H':[[1,1]],'R':[[1]]}])
        self.assertFalse(full_loss_margin(r,F(1,1000000)))
        e=[[F(1)],[-F(1)]]
        self.assertEqual(congruence(r['loss'],e),[[0]])

    def test_rejected_events_contribute_no_loss_and_reset_is_only_a_congruence(self):
        r=word_identity(identity(2),[{'kind':'reset','G':[[1,1],[0,1]]}])
        self.assertEqual(r['loss'],[[0,0],[0,0]])
        self.assertFalse(full_loss_margin(r,F(1,100)))

    def test_nonpositive_noise_cannot_certify_energy(self):
        with self.assertRaises(ValueError):
            word_identity(identity(2),[{'kind':'prediction','F':identity(2),'Q':[[0,1],[1,0]]}])
        with self.assertRaises(ValueError):
            word_identity(identity(2),[{'kind':'correction','H':[[1,0]],'R':[[0]]}])


if __name__=='__main__':
    unittest.main()
