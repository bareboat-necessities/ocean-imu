"""Binary32 tuner tau target/EMA/commit regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T


def cfg():
    return C.CandidateConfig(F(1,10),2,1,1,F(1,10),2,2,1,F(1,100),2,F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def sample(freq=F(1,2)):
    return C.WaveBandSample(freq,True,1,0,False,0,1,1)


class Tests(unittest.TestCase):
    def test_target_horizon_alpha_and_both_ema_evaluation_shapes(self):
        # f=.5 -> tau target=1, sea_time=1, adapt horizon=.4, dt=.01,
        # so rounded x is near .025.  0.9753 is a binary32 exp witness inside
        # the rigorous [1-x,1-x+x^2/2] real enclosure.
        e=B.rn32(F(9753,10000))
        out=T.step(B.rn32(F(1,2)),sample(),cfg(),dt=F(1,100),exp_decay=e)
        self.assertEqual(out.frequency,B.rn32(F(1,2)))
        self.assertEqual(out.tau_target,B.rn32(1))
        self.assertEqual(out.sea_time,B.rn32(1))
        self.assertEqual(out.adapt_sec,B.rn32(F(2,5)))
        self.assertTrue(B.is_binary32(out.alpha))
        self.assertTrue(B.is_binary32(out.next_separate))
        self.assertTrue(B.is_binary32(out.next_fma))
        self.assertEqual(T.committed_tau(out,contracted=False),out.next_separate)
        self.assertEqual(T.committed_tau(out,contracted=True),out.next_fma)

    def test_previous_tau_must_be_actual_stored_binary32(self):
        e=B.rn32(F(9753,10000))
        with self.assertRaisesRegex(ValueError,'actual binary32 stored value'):
            T.step(F(1,3),sample(),cfg(),dt=F(1,100),exp_decay=e)

    def test_detached_exp_witness_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            T.step(B.rn32(F(1,2)),sample(),cfg(),dt=F(1,100),exp_decay=B.rn32(F(9,10)))

    def test_readiness_narrows_remaining_tau_deployment_gaps(self):
        r=T.readiness()
        self.assertTrue(r['shipping_tau_target_EMA_commit_source_shape_matches'])
        self.assertTrue(r['tau_target_float_clamp_graph_materialized'])
        self.assertTrue(r['tau_EMA_separate_mul_add_result_materialized'])
        self.assertTrue(r['tau_EMA_contracted_fma_result_materialized'])
        self.assertTrue(r['pending_commit_passes_stored_tau_directly_to_MEKF_setter'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['shipping_compiler_FP_contraction_mode_qualified'])
        self.assertFalse(r['source_frontend_frequency_binary32_storage_correspondence_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
