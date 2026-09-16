"""Binary32 SeaStateAutoTuner statistics regressions."""
from fractions import Fraction as F
import unittest
from dataclasses import replace

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as X
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E
import test_finite_source_bound_live_word as BASE


def cfg(): return BASE.root_state().runtime.stats_cfg


def coeff(freq=F(1,5),dt=F(1,200)):
    c=cfg(); f=B.rn32(freq); h=B.rn32(dt)
    fe=min(max(f,B.rn32(c.f_min)),B.rn32(c.f_max))
    sea=min(max(B.div(X.HALF,fe),X.TIME_MIN),X.TIME_MAX)
    teff=B.mul(X.TWO,sea)
    req=min(max(B.mul(B.rn32(c.K_periods),teff),B.rn32(c.tau_var_min)),B.rn32(c.tau_var_max))
    lo=X.HORIZON_MIN
    if h>lo: lo=min(h,X.HORIZON_MAX)
    tau=min(max(req,lo),X.HORIZON_MAX)
    arg=B.div(h,tau); elo,ehi=E.exp_minus_enclosure(arg); e=B.rn32((elo+ehi)/2)
    return X.coefficients(c,frequency=f,dt=h,exp_decay=e)


def successor_from_first(env:X.Envelope):
    c=env.coefficients; b=env.before
    return X.State(c.frequency,c.tau_var,env.mean_values[0],env.mean_weights[0],
                   env.sq_values[0],env.sq_weights[0],b.samples+1)


class Tests(unittest.TestCase):
    def test_positive_subnormal_input_is_clamped_before_normal_arithmetic(self):
        floor=coeff(cfg().f_min)
        tiny=F(27,1<<149)
        actual=X.coefficients(cfg(),frequency=tiny,dt=floor.dt,exp_decay=floor.exp_decay)
        self.assertEqual(actual,replace(floor,input_frequency=tiny))
        self.assertEqual(actual.frequency,B.rn32(cfg().f_min))
        with self.assertRaisesRegex(ValueError,'binary32 stats input frequency'):
            X.coefficients(cfg(),frequency=F(1,1<<150),dt=floor.dt,exp_decay=floor.exp_decay)

    def test_first_sample_uses_same_alpha_for_mean_square_and_weights(self):
        s=X.State(); c=coeff(); a=B.rn32(F(1,2)); env=X.envelope(s,c,accel=a); nxt=successor_from_first(env)
        out,cert=X.step(s,c,accel=a,successor=nxt)
        self.assertEqual(out,nxt); self.assertIs(cert.before,s)
        self.assertEqual(out.frequency,c.frequency); self.assertEqual(out.tau_var,c.tau_var)
        self.assertEqual(env.accel_sq,B.mul(a,a))
        self.assertEqual(env.decay,B.sub(B.rn32(1),c.alpha))
        # Shipping DebiasedEMA::isReady() is weight > 1e-6f.  At the default
        # 5 ms update the first alpha already exceeds that threshold, so the
        # first accepted sample is ready even though its debiased variance is 0.
        self.assertTrue(out.var_ready); self.assertEqual(X.variance(out),0)

    def test_successive_sample_must_use_carried_machine_moments(self):
        c=coeff(); a=B.rn32(F(2,5)); s=X.State(); e1=X.envelope(s,c,accel=a); s1=successor_from_first(e1)
        e2=X.envelope(s1,c,accel=a); s2=successor_from_first(e2)
        out,_=X.step(s1,c,accel=a,successor=s2)
        self.assertEqual(out.samples,2); self.assertEqual(e2.before,s1)
        self.assertIn(out.mean_value,e2.mean_values); self.assertIn(out.sq_value,e2.sq_values)

    def test_ready_variance_is_literal_binary32_debiased_readout(self):
        mw=B.rn32(F(1,2)); sw=B.rn32(F(1,2)); mv=B.rn32(F(1,10)); sv=B.rn32(F(1,5))
        s=X.State(B.rn32(F(1,5)),B.rn32(10),mv,mw,sv,sw,100)
        self.assertTrue(s.var_ready)
        mu=B.div(mv,mw); second=B.div(sv,sw)
        expected=max(B.rn32(0),B.sub(second,B.mul(mu,mu)))
        self.assertEqual(X.variance(s),expected)

    def test_detached_successor_fails_closed(self):
        s=X.State(); c=coeff(); a=B.rn32(F(1,2)); env=X.envelope(s,c,accel=a); good=successor_from_first(env)
        bad=X.State(good.frequency,good.tau_var,B.rn32(7),good.mean_weight,good.sq_value,good.sq_weight,good.samples)
        with self.assertRaisesRegex(ValueError,'outside same-step contraction set'):
            X.step(s,c,accel=a,successor=bad)

    def test_detached_exp_witness_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'stats exp witness detached'):
            X.coefficients(cfg(),frequency=B.rn32(F(1,5)),dt=B.rn32(F(1,200)),exp_decay=B.rn32(F(1,2)))

    def test_readiness_closes_machine_stats_topology_not_platform_correspondence(self):
        r=X.readiness()
        for k in ('shipping_tuner_stats_source_shape_matches','frequency_horizon_and_alpha_binary32_graph_materialized',
                  'stats_exp_result_bound_to_same_dt_tau_argument_by_error_profile',
                  'actual_machine_mean_and_square_EMA_successors_bound_to_local_contraction_sets',
                  'mean_and_square_weights_use_same_alpha','binary32_debiased_variance_readout_materialized',
                  'no_exact_real_accel_variance_substituted'):
            self.assertTrue(r[k])
        for k in ('target_exp_libm_correspondence_closed','upstream_machine_band_output_attached',
                  'target_compiler_contraction_membership_closed','startup_frontend_machine_history_attached',
                  'Live_600_step_machine_history_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
