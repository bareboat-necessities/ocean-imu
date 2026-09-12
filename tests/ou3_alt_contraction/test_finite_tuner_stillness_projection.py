"""Exact tuner-relevant StillnessAdapter projection regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_stillness_runtime as FULL
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as T


class Tests(unittest.TestCase):
    def test_moving_projection_ignores_tracker_frequency(self):
        cfg=FULL.Config(energy_alpha=1,energy_thresh=F(1,100))
        s=FULL.State()
        # a/g=1 -> energy=1, so moving regardless of tracker frequency.
        f1=FULL.step(s,cfg,a_vert_up_lp=cfg.gravity,dt=F(1,10),tracker_frequency=F(1,5))
        f2=FULL.step(s,cfg,a_vert_up_lp=cfg.gravity,dt=F(1,10),tracker_frequency=F(4,5))
        p=T.step(T.project(s),cfg,a_vert_up_lp=cfg.gravity,dt=F(1,10))
        self.assertTrue(T.assert_matches_full(p,f1)); self.assertTrue(T.assert_matches_full(p,f2))
        self.assertNotEqual(f1.tracker_frequency_out,f2.tracker_frequency_out)
        self.assertEqual(T.project(f1.state),T.project(f2.state))

    def test_still_prehold_projection_ignores_tracker_frequency(self):
        cfg=FULL.Config()
        s=FULL.State()
        att=FULL.AttenuationWitness(F(9,10))
        f1=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(1,5),attenuation=att)
        f2=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(4,5),attenuation=att)
        p=T.step(T.project(s),cfg,a_vert_up_lp=0,dt=F(1,10),attenuation=att)
        self.assertTrue(T.assert_matches_full(p,f1)); self.assertTrue(T.assert_matches_full(p,f2))
        self.assertEqual(T.project(f1.state),T.project(f2.state))

    def test_posthold_tracker_relaxation_disappears_from_tuner_projection(self):
        cfg=FULL.Config()
        s=FULL.State(energy_ema=0,still_time=F(21,10),freq_init=True,freq_state=F(4,5),last_is_still=True)
        att=FULL.AttenuationWitness(F(1,8))
        # Different tracker values and different frequency relaxation witnesses
        # change only freq_state/output, not energy/still_time/attenuation.
        f1=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(1,5),relax=FULL.RelaxWitness(F(9,10)),attenuation=att)
        f2=FULL.step(s,cfg,a_vert_up_lp=0,dt=F(1,10),tracker_frequency=F(9,10),relax=FULL.RelaxWitness(F(1,2)),attenuation=att)
        p=T.step(T.project(s),cfg,a_vert_up_lp=0,dt=F(1,10),attenuation=att)
        self.assertTrue(T.assert_matches_full(p,f1)); self.assertTrue(T.assert_matches_full(p,f2))
        self.assertNotEqual(f1.tracker_frequency_out,f2.tracker_frequency_out)
        self.assertEqual(T.project(f1.state),T.project(f2.state))

    def test_projection_cannot_consume_attenuation_on_moving_branch(self):
        cfg=FULL.Config(energy_alpha=1,energy_thresh=F(1,100))
        with self.assertRaises(ValueError):
            T.step(T.State(),cfg,a_vert_up_lp=cfg.gravity,dt=F(1,10),attenuation=FULL.AttenuationWitness(1))

    def test_readiness_remains_fail_closed(self):
        r=T.readiness()
        self.assertTrue(r['energy_still_time_projection_independent_of_tracker_frequency'])
        self.assertTrue(r['tracker_algorithm_not_required_for_OU_tuner_stability_word'])
        self.assertFalse(r['stillness_attenuation_exp_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
