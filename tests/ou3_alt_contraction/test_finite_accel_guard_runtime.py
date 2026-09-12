"""Finite AccelVibrationGuard recurrence regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G


class Tests(unittest.TestCase):
    def test_disabled_guard_is_exact_identity_and_consumes_no_witnesses(self):
        s=G.State(); out=G.step(s,G.Config(cutoff_hz=0),(1,2,3),F(1,200))
        self.assertIs(out.state,s); self.assertEqual(out.output,(1,2,3)); self.assertEqual(out.excess_rms,0)
        with self.assertRaises(ValueError):
            G.step(s,G.Config(cutoff_hz=0),(1,2,3),F(1,200),decay=G.DecayWitness(1,1,1,1),rms=G.RmsWitness(0))

    def test_first_enabled_sample_seeds_all_filters_and_is_bit_identity(self):
        out=G.step(G.State(),G.Config(),(1,2,3),F(1,200))
        self.assertTrue(out.state.initialized); self.assertEqual(out.output,(1,2,3)); self.assertEqual(out.state.weight,0)
        self.assertEqual(out.state.stages,((1,2,3),)*4)
        self.assertEqual(out.state.detect_stages,((1,2,3),(0,0,0)))
        self.assertEqual(out.state.removed_ms,(0,0,0))

    def test_two_pole_lowpass_detector_rms_and_partial_blend_are_same_recurrence(self):
        s=G.State(stages=((0,0,0),)*4,detect_stages=((0,0,0),(0,0,0)),removed_ms=(0,0,0),weight=F(1,4),initialized=True)
        # alpha=gamma=1/2, beta=1, slew=1/2, input only on x.
        # low pass: stage0=1, stage1=1/2. detector hp after two stages = 1/2.
        # removed sum=1/4 -> rms=1/2. target clamps to 1; raw weight=5/8.
        out=G.step(s,G.Config(),(2,0,0),F(1,200),decay=G.DecayWitness(F(1,2),F(1,2),1,F(1,2)),rms=G.RmsWitness(F(1,2)))
        self.assertEqual(out.low_passed,(F(1,2),0,0)); self.assertEqual(out.detector_high_passed,(F(1,2),0,0))
        self.assertEqual(out.state.removed_ms,(F(1,4),0,0)); self.assertEqual(out.removed_rms,F(1,2))
        self.assertEqual(out.target,1); self.assertEqual(out.raw_weight,F(5,8)); self.assertEqual(out.state.weight,F(5,8))
        self.assertEqual(out.output,(F(17,16),0,0)); self.assertEqual(out.excess_rms,F(47,100))

    def test_weight_rails_snap_exactly_to_raw_or_lowpass(self):
        s=G.State(stages=((0,0,0),)*4,detect_stages=((0,0,0),(0,0,0)),removed_ms=(0,0,0),weight=F(1,20000),initialized=True)
        # Zero detector target and no slew away -> below epsilon parks at zero.
        out=G.step(s,G.Config(),(0,0,0),F(1,200),decay=G.DecayWitness(1,1,0,0),rms=G.RmsWitness(0))
        self.assertEqual(out.state.weight,0); self.assertEqual(out.output,(0,0,0))
        s2=G.State(stages=((0,0,0),)*4,detect_stages=((0,0,0),(0,0,0)),removed_ms=(1,0,0),weight=1-F(1,20000),initialized=True)
        out2=G.step(s2,G.Config(),(2,0,0),F(1,200),decay=G.DecayWitness(F(1,2),1,0,0),rms=G.RmsWitness(1))
        self.assertEqual(out2.state.weight,1); self.assertEqual(out2.output,out2.low_passed)

    def test_rms_witness_must_belong_to_same_detector_state(self):
        s=G.State(initialized=True)
        with self.assertRaisesRegex(ValueError,'detached'):
            G.step(s,G.Config(),(2,0,0),F(1,200),decay=G.DecayWitness(F(1,2),F(1,2),1,F(1,2)),rms=G.RmsWitness(1))

    def test_readiness_fail_closed_at_binary32_and_complete_word(self):
        r=G.readiness()
        self.assertTrue(r['guard_conditioned_accel_and_excess_RMS_materialized'])
        self.assertTrue(r['guard_engagement_slew_and_rail_snapping_materialized'])
        self.assertFalse(r['guard_exp_sqrt_binary32_ancestry_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
