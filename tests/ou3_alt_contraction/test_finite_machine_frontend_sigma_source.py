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


class Tests(unittest.TestCase):
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
