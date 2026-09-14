"""Persistent actual-machine adaptive-band ledger regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as C
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as Q
from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as X
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E
import test_finite_source_bound_live_word as BASE


def expw(x):
    lo,hi=E.exp_minus_enclosure(x); return B.rn32((lo+hi)/2)

def coeff():
    cfg=BASE.root_state().runtime.band_cfg; f=B.rn32(F(1,5)); h=B.rn32(F(1,200))
    guard=B.div(Q.NYQUIST_FACTOR,h); upper=min(B.rn32(cfg.max_hz),guard)
    low=max(B.rn32(cfg.min_hz),B.mul(B.rn32(cfg.low_ratio),f)); low=min(low,B.div(upper,Q.SPACING))
    high=min(upper,B.mul(B.rn32(cfg.high_ratio),f)); high=max(high,B.mul(low,Q.SPACING)); high=min(high,upper)
    xl=B.mul(B.mul(Q.TWO_PI,low),h); xh=B.mul(B.mul(Q.TWO_PI,high),h)
    return Q.produce(cfg,f_ref=f,dt=h,exp_low=expw(xl),exp_high=expw(xh))

def actual_successor(state,x,c):
    env=C.step(state.machine,x=x,q_low=c.q_low,q_high=c.q_high)
    return C.State(env.lowpass_values[0],env.band_values[0],env.p00_values[0],env.p01_values[0],env.p11_values[0],True)


class Tests(unittest.TestCase):
    def test_active_successor_is_bound_to_same_coefficient_contraction_set(self):
        s=X.initial(); c=coeff(); x=B.rn32(F(1,2)); succ=actual_successor(s,x,c)
        out=X.step(s,c,x=x,successor=succ)
        self.assertEqual(out.state.machine,succ); self.assertEqual(out.state.samples,1)
        self.assertTrue(out.envelope.accepts(succ)); self.assertEqual(out.envelope.q_low,c.q_low)

    def test_successor_outside_set_is_rejected(self):
        s=X.initial(); c=coeff(); x=B.rn32(F(1,2)); bad=C.State(B.rn32(7),B.rn32(0),B.rn32(0),B.rn32(0),B.rn32(0),True)
        with self.assertRaisesRegex(ValueError,'outside same-step'):
            X.step(s,c,x=x,successor=bad)

    def test_readiness_keeps_compiler_membership_and_master_open(self):
        r=X.readiness()
        self.assertTrue(r['actual_machine_band_state_persists_without_contraction_history_branch_explosion'])
        self.assertTrue(r['each_active_successor_must_belong_to_same_step_finite_contraction_set'])
        for k in ('target_compiler_execution_membership_in_local_contraction_relation_closed','target_exp_libm_correspondence_closed',
                  'band_noise_floor_machine_consumer_attached','startup_frontend_machine_history_attached',
                  'Live_600_step_machine_history_attached','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
