"""Common tau/sigma EMA-alpha qualification regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as X
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T
import test_finite_complete_word_tau_qualification as QBASE


def objects(freq=F(1,5),previous=F(11,10)):
    c=QBASE.shipping_runtime().candidate_cfg; d=D.shipping_defaults(qeff_pow_result=B.rn32(1)); dt=B.rn32(F(1,200)); f=B.rn32(freq)
    _,_,_,adapt=T._floats_from_frequency(f,c,dt); x=B.div(dt,adapt); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
    step=T._step_from_frequency(B.rn32(previous),f,c,dt=dt,exp_decay=e)
    return c,d,dt,step


class Tests(unittest.TestCase):
    def test_rederives_tau_step_and_exposes_common_alpha(self):
        c,d,dt,step=objects(); q=X.qualify(step,c,d,dt=dt)
        self.assertIs(q.step,step); self.assertIs(q.candidate_cfg,c); self.assertIs(q.deployment_cfg,d)
        self.assertEqual(q.alpha,step.alpha); self.assertEqual(q.dt,dt)

    def test_detached_common_alpha_scalar_is_rejected(self):
        c,d,dt,step=objects(); bad=replace(c,adapt_tau_sec=c.adapt_tau_sec+1)
        with self.assertRaisesRegex(ValueError,'common-alpha scalars'):
            X.qualify(step,bad,d,dt=dt)

    def test_sigma_target_scale_is_not_misclassified_as_alpha_input(self):
        c,d,dt,step=objects(); changed=replace(c,sigma_coeff=c.sigma_coeff+1,max_sigma=c.max_sigma+1)
        q=X.qualify(step,changed,d,dt=dt)
        self.assertEqual(q.alpha,step.alpha)
        self.assertTrue(X.readiness()['sigma_target_scale_and_clamp_deliberately_deferred_to_sigma_join'])

    def test_constructed_tau_step_with_detached_horizon_is_rejected(self):
        c,d,dt,step=objects()
        bad=replace(step,adapt_sec=B.rn32(F(9,10)))
        with self.assertRaisesRegex(ValueError,'frequency/config horizon'):
            X.qualify(bad,c,d,dt=dt)

    def test_detached_successor_shape_is_rejected(self):
        c,d,dt,step=objects(); changed=B.rn32(step.next_separate+F(1,1<<20))
        self.assertNotEqual(changed,step.next_separate)
        bad=replace(step,next_separate=changed)
        with self.assertRaisesRegex(ValueError,'separate successor'):
            X.qualify(bad,c,d,dt=dt)

    def test_readiness_closes_ancestry_not_libm_or_sigma_ledger(self):
        r=X.readiness()
        for k in ('tau_step_rederived_from_declared_frequency_config_and_dt',
                  'candidate_and_deployment_configs_share_common_tau_sigma_adaptation_scalars',
                  'same_exp_decay_and_alpha_drive_verified_tau_successors',
                  'qualified_alpha_available_for_sigma_machine_recurrence'):
            self.assertTrue(r[k])
        for k in ('target_libm_exp_correspondence_closed','sigma_machine_ledger_attached',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
