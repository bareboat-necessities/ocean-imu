"""Binary32 tuner tau target/EMA/commit regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T


def cfg():
    return C.CandidateConfig(F(1,10),2,1,1,F(1,10),2,2,1,F(1,100),2,F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def sample(freq=F(1,2)):
    return C.WaveBandSample(freq,True,1,0,False,0,1,1)


def exact_candidate(previous,decay):
    return C.step(C.TuneState(previous,F(1,2),F(1,2)),sample(),cfg(),
                  dt=F(1,100),time=F(1,5),last_adapt_time=0,
                  spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(decay,1))


class Tests(unittest.TestCase):
    def test_target_horizon_alpha_and_both_ema_evaluation_shapes(self):
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

    def test_exact_real_shadow_difference_is_retained_as_roundoff_supply(self):
        e=B.rn32(F(9753,10000)); previous=B.rn32(F(1,2))
        binary=T.step(previous,sample(),cfg(),dt=F(1,100),exp_decay=e)
        exact=exact_candidate(previous,e)
        supply=T.roundoff_supply(binary,exact)
        self.assertEqual(supply.residual_separate,
                         binary.next_separate-exact.tune_next.tau_applied)
        self.assertEqual(supply.residual_fma,
                         binary.next_fma-exact.tune_next.tau_applied)

    def test_roundoff_bridge_rejects_detached_exact_frequency(self):
        e=B.rn32(F(9753,10000)); previous=B.rn32(F(1,2))
        binary=T.step(previous,sample(),cfg(),dt=F(1,100),exp_decay=e)
        exact=replace(exact_candidate(previous,e),frequency=F(2,5))
        with self.assertRaisesRegex(ValueError,'frequency detached'):
            T.roundoff_supply(binary,exact)

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
        self.assertTrue(r['source_frontend_frequency_binary32_storage_correspondence_closed'])
        self.assertTrue(r['local_tau_binary32_minus_exact_shadow_supply_exposed'])
        self.assertFalse(r['source_uniform_tau_roundoff_supply_bound_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['shipping_compiler_FP_contraction_mode_qualified'])
        self.assertFalse(r['upstream_WPE_binary32_frequency_production_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
