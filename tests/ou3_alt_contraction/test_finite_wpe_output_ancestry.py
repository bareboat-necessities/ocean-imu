"""Canonical WPE output ancestry regressions; libm correspondence stays open."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_wpe_output_ancestry as X


def cfg(): return W.WPEConfig(1,4,F(1,2),1,180)
def moment_state(log_period=None):
    return W.WPEState(weight=1,elapsed=4,velocity_mean=0,velocity_sq=5,
                      elevation_mean=0,elevation_sq=1,log_period=log_period,
                      last_moment_horizon=4)


class Tests(unittest.TestCase):
    def test_first_output_is_bound_to_new_log_state(self):
        out=X.update(moment_state(),cfg(),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(2,3,7),
                     post=X.CanonicalOutput(7,3,F(1,3)))
        self.assertEqual(out.state.log_period,7)
        self.assertEqual((out.period,out.frequency),(3,F(1,3)))

    def test_existing_current_output_must_name_carried_log_state(self):
        s=moment_state(log_period=2)
        with self.assertRaisesRegex(ValueError,'current canonical output detached'):
            X.update(s,cfg(),current=X.CanonicalOutput(9,5,F(1,5)),
                     post=X.CanonicalOutput(3,1,1),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(2,3,4),
                     log_witness=W.LogUpdateWitness(5,F(1,2),3))

    def test_post_output_must_name_computed_log_successor(self):
        s=moment_state(log_period=2)
        with self.assertRaisesRegex(ValueError,'post canonical output detached'):
            X.update(s,cfg(),current=X.CanonicalOutput(2,5,F(1,5)),
                     post=X.CanonicalOutput(4,1,1),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(2,3,4),
                     log_witness=W.LogUpdateWitness(5,F(1,2),3))

    def test_readiness_closes_ancestry_only(self):
        r=X.readiness()
        self.assertTrue(r['shipping_period_frequency_share_one_log_state_source_shape'])
        self.assertTrue(r['current_output_bound_to_carried_log_period_state'])
        self.assertTrue(r['post_output_bound_to_computed_log_period_successor'])
        self.assertTrue(r['cross_history_reciprocal_output_splicing_forbidden'])
        self.assertFalse(r['period_exp_binary32_libm_correspondence_closed'])
        self.assertFalse(r['frequency_exp_binary32_libm_correspondence_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
