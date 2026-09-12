"""Finite pre-accelerometer Racc runtime regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_racc_runtime as X
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState

ZERO=(F(0),F(0),F(0))

def guard(excess=0):
    s=G.State(initialized=True)
    return G.Result(s,ZERO,ZERO,ZERO,F(excess),F(excess),F(0),F(0))

class Tests(unittest.TestCase):
    def test_deployed_default_dormant_branch_is_literal_no_write(self):
        state=X.State(False,(1,1,1)); cfg=X.Config(vibration_gain=0)
        out=X.step(state,cfg,guard(),nominal_std=(1,1,1),tune=TuneState(1,1,1),
                   preupdate_frequency=F(1,5),live=True)
        self.assertIs(out.state,state); self.assertFalse(out.wrote_Racc)
        self.assertEqual(out.effective_std,(1,1,1)); self.assertEqual(out.scales,(1,1,1))
        self.assertEqual(out.covariance,((1,0,0),(0,1,0),(0,0,1)))

    def test_vibration_only_inflation_uses_same_guard_excess(self):
        state=X.State(False,(3,3,3)); cfg=X.Config(vibration_gain=1)
        out=X.step(state,cfg,guard(4),nominal_std=(3,3,3),tune=TuneState(1,1,1),
                   preupdate_frequency=F(1,5),live=True,
                   effective_sqrt=X.EffectiveSqrtWitness((5,5,5)))
        self.assertTrue(out.wrote_Racc); self.assertTrue(out.state.inflated)
        self.assertEqual(out.excess_rms,4); self.assertEqual(out.effective_std,(5,5,5))
        self.assertEqual(out.covariance,((25,0,0),(0,25,0),(0,0,25)))

    def test_dormant_after_inflation_restores_nominal_exactly_once(self):
        state=X.State(True,(5,5,5)); cfg=X.Config(vibration_gain=1)
        first=X.step(state,cfg,guard(),nominal_std=(3,3,3),tune=TuneState(1,1,1),
                     preupdate_frequency=F(1,5),live=True)
        self.assertTrue(first.wrote_Racc); self.assertTrue(first.restored_nominal)
        self.assertFalse(first.state.inflated); self.assertEqual(first.effective_std,(3,3,3))
        second=X.step(first.state,cfg,guard(),nominal_std=(3,3,3),tune=TuneState(1,1,1),
                      preupdate_frequency=F(1,5),live=True)
        self.assertFalse(second.wrote_Racc); self.assertFalse(second.restored_nominal)

    def test_live_RAO_branch_uses_same_previous_tune_and_frequency(self):
        rao=X.RaoConfig(max_std_scale=2,transition_snr=2,
                        horizontal_x_tau=0,horizontal_y_tau=0,
                        heave_period=1,heave_damping=F(1,2),two_pi=1)
        cfg=X.Config(vibration_gain=0,sigma_coeff=1,rao=rao)
        # sigma_applied=0 -> SNR=0 -> smoothstep weight 1 -> scale=max=2.
        out=X.step(X.State(False,(1,1,1)),cfg,guard(),nominal_std=(1,1,1),
                   tune=TuneState(1,0,1),preupdate_frequency=1,live=True,
                   rao_witness=X.RaoWitness(1,2,2),
                   effective_sqrt=X.EffectiveSqrtWitness((2,2,1)))
        self.assertEqual(out.scales,(2,2,1)); self.assertEqual(out.effective_std,(2,2,1))
        self.assertEqual(out.covariance,((4,0,0),(0,4,0),(0,0,1)))

    def test_RAO_never_applies_before_Live(self):
        rao=X.RaoConfig(max_std_scale=2)
        cfg=X.Config(vibration_gain=0,rao=rao)
        out=X.step(X.State(False,(1,1,1)),cfg,guard(),nominal_std=(1,1,1),
                   tune=TuneState(1,0,1),preupdate_frequency=1,live=False)
        self.assertEqual(out.scales,(1,1,1)); self.assertFalse(out.wrote_Racc)
        with self.assertRaises(ValueError):
            X.step(X.State(False,(1,1,1)),cfg,guard(),nominal_std=(1,1,1),
                   tune=TuneState(1,0,1),preupdate_frequency=1,live=False,
                   rao_witness=X.RaoWitness(1,2,2))

    def test_detached_witnesses_are_rejected(self):
        cfg=X.Config(vibration_gain=1)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.step(X.State(False,(3,3,3)),cfg,guard(4),nominal_std=(3,3,3),
                   tune=TuneState(1,1,1),preupdate_frequency=1,live=True,
                   effective_sqrt=X.EffectiveSqrtWitness((4,5,5)))

    def test_readiness_stays_fail_closed_at_precision_and_stage_ancestry(self):
        r=X.readiness()
        self.assertTrue(r['same_guard_excess_RMS_drives_Racc'])
        self.assertTrue(r['effective_std_to_diagonal_covariance_matches_set_Racc_std'])
        self.assertFalse(r['hypot_and_sqrt_binary32_ancestry_attached'])
        self.assertFalse(r['nominal_Racc_stage_ancestry_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
