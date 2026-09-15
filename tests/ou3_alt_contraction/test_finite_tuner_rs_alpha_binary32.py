"""Binary32 R_S smoothing-horizon / alpha regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as X


# rn32 returns an immutable Fraction; bind the witness defaults once (B008).
MULT_DEFAULT=B.rn32(F(3,2)); TAU_DEFAULT=B.rn32(F(5,2)); DT_DEFAULT=B.rn32(F(1,200))


def exp_witness(mult=MULT_DEFAULT,tau=TAU_DEFAULT,dt=DT_DEFAULT):
    safe=min(max(tau,X.TIME_MIN),X.TIME_MAX)
    requested=B.mul(mult,safe)
    lo=min(max(dt,X.HORIZON_MIN),X.HORIZON_MAX)
    rssec=min(max(requested,lo),X.HORIZON_MAX)
    x=B.div(dt,rssec)
    a,b,_,_=EXP.enclosure(x)
    return B.rn32((a+b)/2)


class Tests(unittest.TestCase):
    def test_deployed_slew_zero_path_materializes_horizon_and_alpha(self):
        mult=B.rn32(F(3,2)); tau=B.rn32(F(5,2)); dt=B.rn32(F(1,200))
        out=X.step(mult=mult,tau_target=tau,dt=dt,exp_decay=exp_witness(mult,tau,dt))
        self.assertEqual(out.safe_tau,tau)
        self.assertEqual(out.requested_horizon,B.mul(mult,tau))
        self.assertEqual(out.RS_sec,out.requested_horizon)
        self.assertEqual(out.exp_argument,B.div(dt,out.RS_sec))
        self.assertEqual(out.alpha,B.sub(B.rn32(1),out.exp_decay))

    def test_tau_time_scale_and_final_horizon_clamps_are_literal(self):
        dt=B.rn32(F(1,200)); mult=B.rn32(F(3,2))
        low=B.rn32(F(1,50)); high=B.rn32(12)
        a=X.step(mult=mult,tau_target=low,dt=dt,exp_decay=exp_witness(mult,low,dt))
        b=X.step(mult=mult,tau_target=high,dt=dt,exp_decay=exp_witness(mult,high,dt))
        self.assertEqual(a.safe_tau,X.TIME_MIN)
        self.assertEqual(b.safe_tau,X.TIME_MAX)
        self.assertEqual(a.RS_sec,B.mul(mult,X.TIME_MIN))
        self.assertEqual(b.RS_sec,B.mul(mult,X.TIME_MAX))

    def test_nonzero_slew_branch_cannot_enter_default_theorem(self):
        with self.assertRaisesRegex(ValueError,'literal slew_log=0'):
            X.step(mult=B.rn32(F(3,2)),tau_target=B.rn32(1),dt=B.rn32(F(1,200)),
                   exp_decay=exp_witness(tau=B.rn32(1)),slew_log=B.rn32(F(1,10)))

    def test_detached_exp_witness_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'SAME rounded'):
            X.step(mult=B.rn32(F(3,2)),tau_target=B.rn32(F(5,2)),dt=B.rn32(F(1,200)),
                   exp_decay=B.rn32(F(1,2)))

    def test_readiness_closes_graph_not_libm_or_master_word(self):
        r=X.readiness()
        for k in ('shipping_RS_horizon_alpha_source_shape_matches',
                  'deployed_slew_log_zero_removes_log_ratio_branch',
                  'dynamic_tau_time_scale_binary32_clamp_materialized',
                  'RS_horizon_binary32_multiply_and_final_clamp_materialized',
                  'rounded_dt_over_RS_horizon_argument_materialized',
                  'RS_exp_witness_bound_to_same_rounded_argument',
                  'one_minus_exp_binary32_alpha_materialized'):
            self.assertTrue(r[k])
        for k in ('target_libm_exp_correspondence_closed','source_uniform_alpha_RS_supply_bound_closed',
                  'persistent_RS_machine_real_recurrence_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
