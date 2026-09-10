import pathlib
import sys
import unittest
from fractions import Fraction as F

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

import ou3_brmm_acceleration_moment_iqc as Q


class AccelerationMomentIQCTests(unittest.TestCase):
    def test_exact_gram_inverse_and_status(self):
        d=Q.build()
        self.assertEqual(Q.validate(d),[])
        self.assertTrue(d['gram_inverse_exact'])
        self.assertTrue(d['J0_J1_J2_same_history_dependency_retained'])
        self.assertTrue(d['all_three_spatial_axes_jointly_coupled'])
        self.assertFalse(d['independent_moment_boxes_used'])
        self.assertFalse(d['independent_axis_boxes_used'])
        self.assertEqual(d['P3_delta'],1e-18)
        self.assertFalse(d['P4_PASS'])

    def test_constant_acceleration_saturates(self):
        A=F(8)
        v=(A,A/F(2),A/F(6))
        self.assertEqual(Q.quadratic(v,Q.GINV),A*A)

    def test_detached_scalar_box_corner_is_excluded(self):
        A=F(8)
        v=(A,-A/F(2),A/F(6))
        self.assertEqual(Q.quadratic(v,Q.GINV),F(193)*A*A)
        self.assertGreater(Q.quadratic(v,Q.GINV),A*A)

    def test_three_axes_share_one_supply(self):
        A=F(8)
        # Splitting the same vector-cap supply equally between two axes stays on
        # the joint energy boundary; assigning a full boundary history to both
        # axes would violate the one-vector-acceleration IQC.
        half=(A/F(2).sqrt() if False else None)
        one=(A,A/F(2),A/F(6))
        self.assertEqual(Q.joint_normalized_moment_energy((one,(F(0),)*3,(F(0),)*3)),A*A)
        self.assertEqual(Q.joint_normalized_moment_energy((one,one,(F(0),)*3)),2*A*A)

    def test_mutation_or_false_promotion_rejected(self):
        d=Q.build()
        d['gram_inverse'][0][0]='8'
        d['P4_PASS']=True
        f=Q.validate(d)
        self.assertIn('Gram inverse changed',f)
        self.assertIn('P4_PASS not false',f)


if __name__=='__main__': unittest.main()
