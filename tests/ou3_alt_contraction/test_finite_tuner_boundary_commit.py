"""Persisted tuner candidate -> next-boundary active parameter regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_tuner_boundary_commit as X
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as PREFIX
from tools.stability.ou3_alt_contraction import finite_tuner_commit as T
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as STILL


def cfg():
    return T.CommitConfig(1,F(1,100),2,F(1,10),False,False,F(1,10),2,1,1,1)


def state(*,pending=True,ready=True,p11=F(1,4)):
    band=B.BandState(p11=p11,ready=ready)
    return PREFIX.State(V.State(),W.WPEState(),band,B.StatsState(),FRONT.LPFState(),
                        STILL.State(),T.TuneState(1,1,1),0,pending,7,F(7,200))


class Tests(unittest.TestCase):
    def test_pending_candidate_uses_carried_band_noise_and_becomes_active(self):
        out=X.apply(state(),cfg(),live=True,bench_noise_sigma=F(3,25),noise_sqrt=B.NoiseSqrtWitness(F(1,2)))
        self.assertFalse(out.state.pending)
        self.assertEqual(out.band_noise_floor_sigma,F(3,50))
        self.assertEqual(out.commit.tau,1); self.assertEqual(out.commit.pseudo_period,F(1,10))
        self.assertEqual(out.commit.Sigma_aw,((1,0,0),(0,1,0),(0,0,1)))
        self.assertEqual(out.commit.R_S,((1,0,0),(0,1,0),(0,0,1)))
        self.assertEqual(out.active.tau,1); self.assertEqual(out.active.Sigma_aw,out.commit.Sigma_aw)
        self.assertEqual(out.active.pseudo_period,F(1,10)); self.assertEqual(out.active.R_S,out.commit.R_S)

    def test_unready_band_uses_raw_bench_noise_without_sqrt(self):
        out=X.apply(state(ready=False,p11=0),cfg(),live=False,bench_noise_sigma=F(3,25))
        self.assertEqual(out.band_noise_floor_sigma,F(3,25)); self.assertIsNone(out.commit.R_S); self.assertIsNone(out.active.R_S)
        with self.assertRaises(ValueError):
            X.apply(state(ready=False,p11=0),cfg(),live=False,bench_noise_sigma=F(3,25),noise_sqrt=B.NoiseSqrtWitness(0))

    def test_ready_band_rejects_detached_covariance_sqrt(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            X.apply(state(),cfg(),live=True,bench_noise_sigma=F(3,25),noise_sqrt=B.NoiseSqrtWitness(F(1,3)))

    def test_no_pending_boundary_is_exact_identity_and_consumes_no_witness(self):
        s=state(pending=False)
        out=X.apply(s,cfg(),live=True,bench_noise_sigma=F(3,25))
        self.assertIs(out.state,s); self.assertIsNone(out.commit); self.assertIsNone(out.active); self.assertIsNone(out.band_noise_floor_sigma)
        with self.assertRaisesRegex(ValueError,'no-pending'):
            X.apply(s,cfg(),live=True,bench_noise_sigma=F(3,25),noise_sqrt=B.NoiseSqrtWitness(F(1,2)))

    def test_readiness_fail_closed_only_at_roundoff_and_complete_word(self):
        r=X.readiness()
        self.assertTrue(r['candidate_to_active_parameter_ancestry_closed'])
        self.assertTrue(r['band_noise_floor_derived_from_carried_band_state'])
        self.assertFalse(r['commit_sqrt_and_binary32_roundoff_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
