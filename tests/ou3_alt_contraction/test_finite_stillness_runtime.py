"""Tracker stillness finite runtime regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as S


def cfg(): return S.Config(gravity=1,energy_alpha=F(1,20),energy_thresh=F(8,10000),still_thresh=2,relax_tau=1,target_freq=F(1,5))


class Tests(unittest.TestCase):
    def test_moving_branch_resets_hold_and_passes_tracker_frequency(self):
        s=S.State(energy_ema=0,still_time=5,freq_init=True,freq_state=F(1,4),last_is_still=True)
        out=S.step(s,cfg(),a_vert_up_lp=1,dt=F(1,10),tracker_frequency=F(2,5))
        self.assertFalse(out.state.last_is_still)
        self.assertEqual(out.state.still_time,0)
        self.assertEqual(out.tracker_frequency_out,F(2,5))
        self.assertEqual(out.variance_attenuation,1)

    def test_still_pre_hold_tracks_input_without_relaxation(self):
        out=S.step(S.State(),cfg(),a_vert_up_lp=0,dt=1,tracker_frequency=F(2,5),attenuation=S.AttenuationWitness(F(1,2)))
        self.assertTrue(out.state.last_is_still)
        self.assertEqual(out.state.still_time,1)
        self.assertEqual(out.tracker_frequency_out,F(2,5))
        self.assertEqual(out.variance_attenuation,F(1,2))

    def test_still_post_hold_relaxes_and_caps_time(self):
        s=S.State(energy_ema=0,still_time=59,freq_init=True,freq_state=F(2,5),last_is_still=True)
        out=S.step(s,cfg(),a_vert_up_lp=0,dt=3,tracker_frequency=F(1,2),relax=S.RelaxWitness(F(1,2)),attenuation=S.AttenuationWitness(F(1,4)))
        self.assertEqual(out.state.still_time,60)
        self.assertEqual(out.tracker_frequency_out,F(3,10))
        self.assertEqual(out.variance_attenuation,F(1,4))

    def test_branch_witnesses_cannot_be_used_on_wrong_branch(self):
        with self.assertRaises(ValueError):
            S.step(S.State(),cfg(),a_vert_up_lp=1,dt=1,tracker_frequency=F(1,2),attenuation=S.AttenuationWitness(F(1,2)))
        with self.assertRaises(ValueError):
            S.step(S.State(),cfg(),a_vert_up_lp=0,dt=1,tracker_frequency=F(1,2),relax=S.RelaxWitness(F(1,2)),attenuation=S.AttenuationWitness(F(1,2)))

    def test_readiness_remains_fail_closed_upstream(self):
        r=S.readiness(); self.assertTrue(r['energy_ema_and_stillness_threshold_materialized']); self.assertTrue(r['variance_attenuation_branch_materialized']); self.assertFalse(r['tracker_frequency_and_vertical_input_same_history_attached']); self.assertFalse(r['relax_and_attenuation_exp_binary32_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
