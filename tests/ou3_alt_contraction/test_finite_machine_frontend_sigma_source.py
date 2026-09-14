"""Same-sample machine band/statistics sigma-source regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as BR
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as BC
from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as X
from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as ST
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
import test_finite_source_bound_live_word as BASE


def expw(x):
    lo,hi=E.exp_minus_enclosure(x); return B.rn32((lo+hi)/2)

def sqrtw(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)

def band_coeff():
    c=BASE.root_state().runtime.band_cfg; f=B.rn32(F(1,5)); h=B.rn32(F(1,200))
    guard=B.div(BC.NYQUIST_FACTOR,h); upper=min(B.rn32(c.max_hz),guard)
    low=max(B.rn32(c.min_hz),B.mul(B.rn32(c.low_ratio),f)); low=min(low,B.div(upper,BC.SPACING))
    high=min(upper,B.mul(B.rn32(c.high_ratio),f)); high=max(high,B.mul(low,BC.SPACING)); high=min(high,upper)
    xl=B.mul(B.mul(BC.TWO_PI,low),h); xh=B.mul(B.mul(BC.TWO_PI,high),h)
    return BC.produce(c,f_ref=f,dt=h,exp_low=expw(xl),exp_high=expw(xh))
def stats_coeff():
    c=BASE.root_state().runtime.stats_cfg; f=B.rn32(F(1,5)); h=B.rn32(F(1,200))
    fe=min(max(f,B.rn32(c.f_min)),B.rn32(c.f_max)); sea=min(max(B.div(ST.HALF,fe),ST.TIME_MIN),ST.TIME_MAX)
    te=B.mul(ST.TWO,sea); req=min(max(B.mul(B.rn32(c.K_periods),te),B.rn32(c.tau_var_min)),B.rn32(c.tau_var_max))
    lo=ST.HORIZON_MIN
    if h>lo: lo=min(h,ST.HORIZON_MAX)
    tau=min(max(req,lo),ST.HORIZON_MAX); arg=B.div(h,tau)
    return ST.coefficients(c,frequency=f,dt=h,exp_decay=expw(arg))

def build_successors(state:X.State,x):
    bc=band_coeff(); env=BR.step(state.band.machine,x=x,q_low=bc.q_low,q_high=bc.q_high)
    band=BR.State(env.lowpass_values[0],env.band_values[0],env.p00_values[0],env.p01_values[0],env.p11_values[0],True)
    sc=stats_coeff(); senv=ST.envelope(state.stats,sc,accel=band.band)
    stats=ST.State(sc.frequency,sc.tau_var,senv.mean_values[0],senv.mean_weights[0],senv.sq_values[0],senv.sq_weights[0],state.stats.samples+1)
    sg=sqrtw(band.p11)
    return bc,band,sc,stats,sg


def source_step(state,*,band_cfg,stats_cfg,frequency,dt=F(1,200),bench_noise_sigma=0,x=0,last=False):
    """Conditional finite-machine regression witness, not native libm evidence."""
    f=B.rn32(frequency); h=B.rn32(dt)
    fr=state.stats.frequency if state.stats.frequency is not None and state.stats.frequency>0 else f
    fr=min(max(fr,B.rn32(band_cfg.tune_freq_floor)),B.rn32(band_cfg.tune_freq_ceil))
    upper=min(B.rn32(band_cfg.max_hz),B.div(BC.NYQUIST_FACTOR,h))
    low=min(max(B.rn32(band_cfg.min_hz),B.mul(B.rn32(band_cfg.low_ratio),fr)),B.div(upper,BC.SPACING))
    high=min(max(min(upper,B.mul(B.rn32(band_cfg.high_ratio),fr)),B.mul(low,BC.SPACING)),upper)
    bc=BC.produce(band_cfg,f_ref=fr,dt=h,
        exp_low=expw(B.mul(B.mul(BC.TWO_PI,low),h)),exp_high=expw(B.mul(B.mul(BC.TWO_PI,high),h)))
    env=BR.step(state.band.machine,x=B.rn32(x),q_low=bc.q_low,q_high=bc.q_high)
    i=-1 if last else 0
    bn=BR.State(env.lowpass_values[i],env.band_values[i],env.p00_values[i],env.p01_values[i],env.p11_values[i],True)
    fe=min(max(f,B.rn32(stats_cfg.f_min)),B.rn32(stats_cfg.f_max))
    sea=min(max(B.div(ST.HALF,fe),ST.TIME_MIN),ST.TIME_MAX)
    req=min(max(B.mul(B.rn32(stats_cfg.K_periods),B.mul(ST.TWO,sea)),B.rn32(stats_cfg.tau_var_min)),B.rn32(stats_cfg.tau_var_max))
    lo=ST.HORIZON_MIN if h<=ST.HORIZON_MIN else min(h,ST.HORIZON_MAX)
    tau=min(max(req,lo),ST.HORIZON_MAX)
    sc=ST.coefficients(stats_cfg,frequency=f,dt=h,exp_decay=expw(B.div(h,tau)))
    se=ST.envelope(state.stats,sc,accel=bn.band)
    sn=ST.State(sc.frequency,sc.tau_var,se.mean_values[i],se.mean_weights[i],se.sq_values[i],se.sq_weights[i],state.samples+1)
    return X.step(state,band_coefficients=bc,band_input=B.rn32(x),band_successor=bn,
        stats_coefficients=sc,stats_successor=sn,bench_noise_sigma=B.rn32(bench_noise_sigma),
        noise_sqrt_gain=sqrtw(bn.p11),accel_variance=ST.variance_outcomes(sn)[i])


def ready_pair(*,samples=0):
    from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as L
    def one(gain):
        return X.State(L.State(BR.State(p11=F(gain),ready=True),samples),ST.State(samples=samples))
    return X.Pair(one(121),one(144))


class Tests(unittest.TestCase):
    def test_two_sample_source_join_retains_lagged_band_frequency(self):
        from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as W
        from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WP
        from dataclasses import replace
        word=BASE.root_state(); runtime=word.runtime
        c=runtime.band_cfg; sc=runtime.stats_cfg
        first=source_step(X.initial(),band_cfg=c,stats_cfg=sc,frequency=B.rn32(F(1,5)),x=F(1,2))
        previous=first.state
        shadow=WP.WPEState(log_period=F(7,10),usable_period=True)
        getter=W.getters(W.bind_log_state(shadow,B.rn32(shadow.log_period)),period_exp=2,frequency_exp=F(1,2))
        fr=W.through_statistics(W.tuner_frequency(shadow,min_hz=B.rn32(c.tune_freq_floor),max_hz=B.rn32(c.tune_freq_ceil),
             getter=getter,shadow_frequency=F(1,2)),sc,exact_min_hz=c.tune_freq_floor,exact_max_hz=c.tune_freq_ceil)
        out=source_step(previous,band_cfg=c,stats_cfg=sc,frequency=F(1,2),x=F(3,4))
        self.assertIs(X.bind_step(previous,out,frequency=fr,band_cfg=c,stats_cfg=sc,dt=F(1,200),bench_noise_sigma=0),out)
        self.assertEqual(out.band.coefficients.f_ref,previous.stats.frequency)
        self.assertNotEqual(out.band.coefficients.f_ref,fr.stats_stored.stored_hz)
        self.assertEqual(out.state.stats.frequency,fr.stats_stored.stored_hz)
        with self.assertRaisesRegex(ValueError,'history detached'):
            X.bind_step(X.initial(),out,frequency=fr,band_cfg=c,stats_cfg=sc,dt=F(1,200),bench_noise_sigma=0)
        with self.assertRaisesRegex(ValueError,'re-executed source'):
            X.bind_step(previous,out,frequency=fr,band_cfg=c,stats_cfg=sc,dt=F(1,200),bench_noise_sigma=F(1,10))
        with self.assertRaisesRegex(ValueError,'SAME rounded argument'):
            X.bind_step(previous,out,frequency=fr,band_cfg=replace(c,low_ratio=c.low_ratio/F(2)),stats_cfg=sc,dt=F(1,200),bench_noise_sigma=0)

    def test_variance_readout_retains_separate_and_contracted_subtraction(self):
        mu=B.rn32(1-F(3,1<<24))
        state=ST.State(mean_value=mu,mean_weight=1,sq_value=1,sq_weight=1)
        outcomes=ST.variance_outcomes(state)
        self.assertEqual(len(outcomes),2)
        self.assertIn(ST.variance(state),outcomes)
        self.assertIn(B.fma(-mu,mu,1),outcomes)
        self.assertTrue(all(B.is_binary32(x) and x>=0 for x in outcomes))

    def test_unconsumed_boundary_does_not_read_new_or_old_band_covariance(self):
        out=X.boundary_floors(ready_pair(),bench_noise_sigma=F(1,10),required=False)
        self.assertEqual(out,(None,None))
        with self.assertRaisesRegex(ValueError,'unconsumed boundary'):
            X.boundary_floors(ready_pair(),bench_noise_sigma=F(1,10),required=False,separate_sqrt_gain=11)

    def test_same_new_band_output_feeds_stats_variance_and_noise(self):
        s=X.initial(); x=B.rn32(F(1,2)); bc,bn,sc,sn,sg=build_successors(s,x)
        out=X.step(s,band_coefficients=bc,band_input=x,band_successor=bn,
                   stats_coefficients=sc,stats_successor=sn,bench_noise_sigma=B.rn32(F(3,100)),noise_sqrt_gain=sg)
        self.assertEqual(out.stats_envelope.accel,out.state.band.machine.band)
        self.assertEqual(out.accel_variance,ST.variance(out.state.stats))
        self.assertEqual(out.band_noise_sigma,out.noise.noise_sigma)
        self.assertEqual(out.state.samples,1)

    def test_stats_successor_from_different_band_output_is_rejected(self):
        s=X.initial(); x=B.rn32(F(1,2)); bc,bn,sc,sn,sg=build_successors(s,x)
        # Replace stats values with a valid local successor generated from a
        # different accel; it must not splice into the actual band successor.
        wrong_env=ST.envelope(s.stats,sc,accel=B.rn32(F(3,4)))
        wrong=ST.State(sc.frequency,sc.tau_var,wrong_env.mean_values[0],wrong_env.mean_weights[0],wrong_env.sq_values[0],wrong_env.sq_weights[0],1)
        with self.assertRaisesRegex(ValueError,'outside same-step contraction set'):
            X.step(s,band_coefficients=bc,band_input=x,band_successor=bn,
                   stats_coefficients=sc,stats_successor=wrong,bench_noise_sigma=B.rn32(F(3,100)),noise_sqrt_gain=sg)

    def test_readiness_closes_machine_variance_noise_source_not_stillness_or_platform(self):
        r=X.readiness()
        for k in ('persistent_machine_band_and_stats_states_composed','statistics_consumes_same_new_machine_band_output',
                  'machine_accel_variance_derived_from_same_stats_successor','machine_band_noise_floor_derived_from_same_band_successor',
                  'band_and_stats_sample_counts_advance_together','external_stats_frequency_not_falsely_identified_with_band_reference_frequency',
                  'binary32_sigma_target_can_consume_machine_variance_and_noise'):
            self.assertTrue(r[k])
        for k in ('stillness_machine_state_attached','WPE_and_band_stats_libm_correspondence_closed',
                  'target_compiler_contraction_membership_closed','startup_frontend_machine_history_attached',
                  'Live_600_step_machine_history_attached','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
