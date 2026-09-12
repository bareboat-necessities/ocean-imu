"""Staged tuner commit regressions; not frontend/source qualification."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_tuner_commit as T


def cfg(cubic=True,scaled=True):
    return T.CommitConfig(F(3,2),F(1,100),2,F(3,200),scaled,cubic,F(2,5),5,F(4,5),F(6,5),F(7,5))

class Tests(unittest.TestCase):
    def test_pending_false_has_no_commit(self):
        self.assertIsNone(T.commit(T.TuneState(1,F(1,10),1),cfg(),pending=False,live=True,band_noise_floor_sigma=F(1,20),rs_sqrt_scale=F(1,10)))

    def test_same_tau_drives_clamped_cadence(self):
        t=T.TuneState(1,F(1,10),1); c=cfg(False,True); out=T.commit(t,c,pending=True,live=False,band_noise_floor_sigma=F(1,20))
        self.assertEqual(out.tau,1); self.assertEqual(out.pseudo_period,F(3,2))
        low=T.commit(T.TuneState(F(1,1000),F(1,10),1),c,pending=True,live=False,band_noise_floor_sigma=F(1,20)); self.assertEqual(low.pseudo_period,c.pseudo_period_min)
        high=T.commit(T.TuneState(2,F(1,10),1),c,pending=True,live=False,band_noise_floor_sigma=F(1,20)); self.assertEqual(high.pseudo_period,c.pseudo_period_max)

    def test_sigma_floor_and_horizontal_factor_build_one_stationary_covariance(self):
        t=T.TuneState(1,F(1,100),1); out=T.commit(t,cfg(False),pending=True,live=False,band_noise_floor_sigma=F(3,50))
        # floor=max(.05,.06)=.06; sH=.8*.06=.048
        self.assertEqual(out.Sigma_aw[2][2],F(9,2500)); self.assertEqual(out.Sigma_aw[0][0],F(144,62500)); self.assertEqual(out.Sigma_aw[1][1],out.Sigma_aw[0][0])

    def test_live_RS_uses_same_realized_period_and_anisotropic_std_factors(self):
        # period=3/2, fixed=3/200 -> sqrt ratio = 1/10 exactly.
        t=T.TuneState(1,F(1,10),2); c=cfg(True,True); out=T.commit(t,c,pending=True,live=True,band_noise_floor_sigma=F(1,20),rs_sqrt_scale=F(1,10))
        z=F(1,5)  # base 2 * info .1
        self.assertEqual(out.R_S,((z*z*F(36,25),0,0),(0,z*z*F(49,25),0),(0,0,z*z)))

    def test_wrong_information_sqrt_witness_rejected(self):
        with self.assertRaises(ValueError): T.commit(T.TuneState(1,F(1,10),1),cfg(True,True),pending=True,live=True,band_noise_floor_sigma=F(1,20),rs_sqrt_scale=F(1,9))

    def test_nonlive_commit_does_not_apply_RS(self):
        out=T.commit(T.TuneState(1,F(1,10),1),cfg(True,True),pending=True,live=False,band_noise_floor_sigma=F(1,20))
        self.assertIsNone(out.R_S)

    def test_periodic_online_commit_does_not_queue_aw_floor(self):
        t=T.TuneState(1,F(1,10),1); c=cfg(False); a=T.commit(t,c,pending=True,live=False,band_noise_floor_sigma=F(1,20),sync_covariance=False); b=T.commit(t,c,pending=True,live=False,band_noise_floor_sigma=F(1,20),sync_covariance=True)
        self.assertIsNone(a.aw_floor_target); self.assertEqual(b.aw_floor_target,b.Sigma_aw)

    def test_readiness_stays_fail_closed(self):
        r=T.readiness(); self.assertTrue(r['same_tau_drives_OU_and_S_cadence']); self.assertTrue(r['live_RS_from_same_TuneState_and_realized_period']); self.assertFalse(r['WPE_band_sigma_candidate_history_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
