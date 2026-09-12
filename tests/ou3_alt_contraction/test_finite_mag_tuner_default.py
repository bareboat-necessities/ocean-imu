"""Default MagAutoTuner finite state-machine regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as X
from tools.stability.ou3_alt_contraction import finite_mag_accumulator as ACC
from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as GAUGE
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT


def root(v): return X.SqrtWitness(v*v,v)
def troot(v): return TILT.SqrtWitness(v*v,v)


class Tests(unittest.TestCase):
    def test_ready_latch_is_exact_identity_and_consumes_no_new_witnesses(self):
        acc=ACC.State((3,4,12),13,1,F(1,200),1)
        s=X.State(acc,0,True,(3,4,12),(5,0,12),1,(3,4,12))
        out=X.step(s,X.Config(),q_tilt_bw=(0,0,0,0),mag_body=(0,0,0),dt=0)
        self.assertTrue(out.returned_ready); self.assertFalse(out.accepted); self.assertEqual(out.state,s)
        with self.assertRaisesRegex(ValueError,'consumes no new arithmetic'):
            X.step(s,X.Config(),q_tilt_bw=(1,0,0,0),mag_body=(3,4,0),dt=0,q_norm=root(1))

    def test_invalid_quaternion_and_small_mag_reject_without_accumulation(self):
        cfg=X.Config()
        out=X.step(X.State(),cfg,q_tilt_bw=(0,0,0,0),mag_body=(3,4,0),dt=F(1,200),
                   q_norm=root(0),mag_norm=root(5))
        self.assertEqual(out.rejection,'quaternion_norm'); self.assertEqual(out.state.accumulator,ACC.State())
        out=X.step(X.State(),cfg,q_tilt_bw=(1,0,0,0),mag_body=(0,0,F(1,2000)),dt=F(1,200),
                   q_norm=root(1),mag_norm=root(F(1,2000)))
        self.assertEqual(out.rejection,'mag_norm'); self.assertEqual(out.state.rejected_count,1)

    def test_accepted_sample_is_rotated_by_same_normalized_tilt_quaternion(self):
        out=X.step(X.State(),X.Config(min_samples=2,min_window=0),
                   q_tilt_bw=(1,0,0,0),mag_body=(3,4,0),dt=F(1,200),
                   q_norm=root(1),mag_norm=root(5))
        self.assertTrue(out.accepted); self.assertFalse(out.returned_ready)
        self.assertEqual(out.state.accumulator.world_sum,(3,4,0))
        self.assertEqual(out.state.accumulator.norm_sum,5)
        self.assertEqual(out.state.last_world_sample,(3,4,0)); self.assertEqual(out.state.last_sample_weight,1)

    def test_boat_quaternion_entry_removes_yaw_before_same_tuner_step(self):
        # Pure yaw q=(3/5,0,0,4/5) is stripped to identity before accumulation.
        q=(F(3,5),0,0,F(4,5))
        yaw=TILT.YawHalfWitness(F(-7,25),F(24,25),troot(1),F(3,5),F(4,5))
        out=X.step_from_boat_quaternion(X.State(),X.Config(min_samples=2,min_window=0),
                    q_boat_bw=q,mag_body=(3,4,0),dt=F(1,200),
                    boat_q_norm=troot(1),yaw_half=yaw,mag_norm=root(5))
        self.assertTrue(out.accepted); self.assertEqual(out.state.last_world_sample,(3,4,0))
        self.assertEqual(out.state.accumulator.world_sum,(3,4,0))

    def test_running_norm_outlier_is_strict_and_does_not_mutate_accepted_statistics(self):
        a=ACC.State((3,4,0),5,1,F(1,200),1)
        s=X.State(a)
        rejected=X.step(s,X.Config(min_samples=10,min_window=0),
                        q_tilt_bw=(1,0,0,0),mag_body=(7,0,0),dt=F(1,200),
                        q_norm=root(1),mag_norm=root(7))
        self.assertEqual(rejected.rejection,'running_norm_outlier')
        self.assertEqual(rejected.state.accumulator,a)
        boundary=F(27,4)
        accepted=X.step(s,X.Config(min_samples=10,min_window=0),
                        q_tilt_bw=(1,0,0,0),mag_body=(boundary,0,0),dt=F(1,200),
                        q_norm=root(1),mag_norm=root(boundary))
        self.assertTrue(accepted.accepted); self.assertIsNone(accepted.rejection)

    def test_nonpositive_dt_falls_back_to_configured_period(self):
        cfg=X.Config(min_samples=2,min_window=0,sample_dt=F(1,100))
        out=X.step(X.State(),cfg,q_tilt_bw=(1,0,0,0),mag_body=(3,4,0),dt=0,
                   q_norm=root(1),mag_norm=root(5))
        self.assertEqual(out.state.accumulator.accepted_window,F(1,100))

    def test_finalize_builds_same_mean_and_gauge_fixed_reference(self):
        cfg=X.Config(min_samples=1,min_window=0)
        out=X.step(X.State(),cfg,q_tilt_bw=(1,0,0,0),mag_body=(3,4,12),dt=F(1,200),
                   q_norm=root(1),mag_norm=root(13),mean_norm=root(13),
                   horizontal_sqrt=GAUGE.HorizontalSqrt(25,5))
        self.assertTrue(out.accepted); self.assertTrue(out.returned_ready); self.assertTrue(out.state.ready)
        self.assertEqual(out.state.mean,(3,4,12)); self.assertEqual(out.state.world_reference,(5,0,12))

    def test_horizontal_fraction_gate_precedes_ready_latch(self):
        cfg=X.Config(min_samples=1,min_window=0,min_horizontal_fraction=F(7,10))
        out=X.step(X.State(),cfg,q_tilt_bw=(1,0,0,0),mag_body=(3,0,4),dt=F(1,200),
                   q_norm=root(1),mag_norm=root(5),mean_norm=root(5),
                   horizontal_sqrt=GAUGE.HorizontalSqrt(9,3))
        self.assertTrue(out.accepted); self.assertFalse(out.returned_ready)
        self.assertEqual(out.rejection,'horizontal_fraction'); self.assertFalse(out.state.ready)

    def test_optional_weighting_and_hard_iron_cannot_enter_default_proof_branch(self):
        with self.assertRaises(NotImplementedError):
            X.step(X.State(),X.Config(quality_weighting=True),q_tilt_bw=(1,0,0,0),mag_body=(3,4,0),dt=1,
                   q_norm=root(1),mag_norm=root(5))
        with self.assertRaises(NotImplementedError):
            X.step(X.State(),X.Config(estimate_hard_iron=True),q_tilt_bw=(1,0,0,0),mag_body=(3,4,0),dt=1,
                   q_norm=root(1),mag_norm=root(5))

    def test_readiness_keeps_binary32_source_and_optional_branches_open(self):
        r=X.readiness()
        self.assertTrue(r['raw_sample_acceptance_boolean_removed'])
        self.assertTrue(r['startup_tilt_quaternion_runtime_ancestry_attached'])
        self.assertFalse(r['startup_tilt_atan2_AngleAxis_binary32_attached'])
        self.assertFalse(r['startup_mag_sensor_source_ancestry_attached'])
        self.assertFalse(r['optional_hard_iron_branch_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
