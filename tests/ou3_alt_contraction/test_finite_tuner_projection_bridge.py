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


def front():
    return B.FrontendResult(B.BandState(),B.StatsState(frequency=F(1,2)),F(1,2),F(1,2),F(1,2),0,True,4,0)


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

    def test_readiness_fail_closed(self):
        r=BR.readiness()
        self.assertTrue(r['no_tracker_frequency_operand_in_tuner_candidate'])
        self.assertTrue(r['dominant_frequency_tracker_removed_from_OU_stability_interconnection'])
        self.assertFalse(r['sigma_wave_sqrt_transcendental_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
