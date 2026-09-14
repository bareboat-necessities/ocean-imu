"""Tracker-free tuner bridge regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as FULL
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as P
from tools.stability.ou3_alt_contraction import finite_tuner_projection_bridge as BR
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as IC
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def front(freq=F(1,2),variance=F(4)):
    return B.FrontendResult(B.BandState(),B.StatsState(frequency=freq),freq,freq,freq,0,True,variance,0)


def projected_still():
    cfg=FULL.Config(); s=FULL.State(energy_ema=0,still_time=1,freq_init=True,freq_state=F(1,2),last_is_still=True)
    return P.step(P.project(s),cfg,a_vert_up_lp=0,dt=F(1,10),attenuation=FULL.AttenuationWitness(F(1,4)))


def shipping_interval_cfg():
    return C.CandidateConfig(
        TARGET.FLOOR,TARGET.CEIL,F(1),F(1),TARGET.TAU_MIN,TARGET.TAU_MAX,F(4),
        B32.div(B32.rn32(F(3,200)),B32.rn32(F(11,10))),B32.rn32(F(1,200)),B32.rn32(F(3,20)),
        F(1,100),F(100),F(1),F(1,2),F(1),B32.rn32(F(9,5)),B32.rn32(F(2,5)),F(3,2),0,F(1,10),True)


class Tests(unittest.TestCase):
    def test_projected_sample_equals_full_shipping_sample(self):
        cfg=FULL.Config(); s=FULL.State(energy_ema=0,still_time=1,freq_init=True,freq_state=F(1,2),last_is_still=True)
        att=FULL.AttenuationWitness(F(1,4))
        full=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(9,10),attenuation=att)
        proj=P.step(P.project(s),cfg,a_vert_up_lp=0,dt=F(1,10),attenuation=att)
        a=BR.sample_from_projection(front(),proj,sigma_wave_sqrt=1)
        b=C.sample_from_runtime(front(),full,sigma_wave_sqrt=1)
        self.assertEqual(a,b)

    def test_two_different_tracker_paths_have_identical_tuner_sample(self):
        cfg=FULL.Config(); s=FULL.State(energy_ema=0,still_time=F(21,10),freq_init=True,freq_state=F(4,5),last_is_still=True)
        att=FULL.AttenuationWitness(F(1,4))
        full1=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(1,5),relax=FULL.RelaxWitness(F(9,10)),attenuation=att)
        full2=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(9,10),relax=FULL.RelaxWitness(F(1,2)),attenuation=att)
        proj=P.step(P.project(s),cfg,a_vert_up_lp=0,dt=F(1,10),attenuation=att)
        projected=BR.sample_from_projection(front(),proj,sigma_wave_sqrt=1)
        self.assertNotEqual(full1.tracker_frequency_out,full2.tracker_frequency_out)
        self.assertEqual(projected,C.sample_from_runtime(front(),full1,sigma_wave_sqrt=1))
        self.assertEqual(projected,C.sample_from_runtime(front(),full2,sigma_wave_sqrt=1))

    def test_interval_projection_accepts_general_shipping_spectral_cell(self):
        # Use f=.2 and sigma=1 on the same tracker-free projection.  The
        # resulting shipping cadence has irrational sqrt(T_S), so the legacy
        # exact-rational SpectralWitness cannot represent it; the interval path
        # must do so without choosing a fake root.
        still=projected_still(); f=front(F(1,5),F(1))
        prev=IC.IntervalTuneState.point(TuneState(F(11,10),F(1,2),F(1,2)))
        out=BR.step_interval(prev,f,still,shipping_interval_cfg(),sigma_wave_sqrt=F(1,2),
                             dt=F(1,200),time=F(1,5),last_adapt_time=0,
                             ema=C.EmaWitness(F(99,100),F(49,50)))
        self.assertEqual(out.target.frequency,F(1,5))
        self.assertEqual(out.target.tau_target,F(5,2))
        self.assertLessEqual(out.tune_next.RS_lo,out.tune_next.RS_hi)

    def test_readiness_fail_closed(self):
        r=BR.readiness()
        self.assertTrue(r['no_tracker_frequency_operand_in_tuner_candidate'])
        self.assertTrue(r['dominant_frequency_tracker_removed_from_OU_stability_interconnection'])
        self.assertTrue(r['general_SpectralMSE_interval_candidate_projection_available'])
        self.assertTrue(r['actual_0p2Hz_prior_requires_no_fake_rational_spectral_witness'])
        self.assertFalse(r['sigma_wave_sqrt_transcendental_attached'])
        self.assertFalse(r['binary32_sqrt_pow_exp_correspondence_closed'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
