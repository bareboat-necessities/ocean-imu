"""Exact-real interval SpectralMSE commit regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C
from tools.stability.ou3_alt_contraction import finite_tuner_interval_commit as X


def cfg():
    return C.CommitConfig(F(3,2),F(1,100),2,F(3,200),False,False,
                          F(2,5),5,F(4,5),F(6,5),F(7,5))


class Tests(unittest.TestCase):
    def test_live_commit_propagates_RS_interval_without_selecting_point(self):
        t=I.IntervalTuneState(F(11,10),F(1,2),F(1,2),F(3,4))
        out=X.commit(t,cfg(),pending=True,live=True,band_noise_floor_sigma=0)
        self.assertEqual(out.tau,t.tau_applied)
        self.assertEqual(out.pseudo_period,C.cadence(C.TuneState(t.tau_applied,t.sigma_applied,t.RS_lo),cfg()))
        self.assertIsNotNone(out.RS_cov_lo); self.assertIsNotNone(out.RS_cov_hi)
        for i in range(3): self.assertLessEqual(out.RS_cov_lo[i][i],out.RS_cov_hi[i][i])
        self.assertFalse(out.pending_after)

    def test_nonlive_commit_carries_no_RS_covariance(self):
        t=I.IntervalTuneState(F(1),F(1),F(1,4),F(1,2))
        out=X.commit(t,cfg(),pending=True,live=False,band_noise_floor_sigma=0)
        self.assertIsNone(out.RS_cov_lo); self.assertIsNone(out.RS_cov_hi)

    def test_no_pending_is_literal_noop_and_cubic_branch_stays_fail_closed(self):
        t=I.IntervalTuneState(F(1),F(1),F(1,4),F(1,2))
        self.assertIsNone(X.commit(t,cfg(),pending=False,live=True,band_noise_floor_sigma=0))
        c=cfg(); cubic=C.CommitConfig(c.pseudo_tau_ratio,c.pseudo_period_min,c.pseudo_period_max,
             c.pseudo_fixed_period,c.tau_scaled_cadence,True,c.min_R_S,c.max_R_S,
             c.S_factor,c.R_S_x_factor,c.R_S_y_factor)
        with self.assertRaisesRegex(ValueError,'non-Cubic SpectralMSE'):
            X.commit(t,cubic,pending=True,live=True,band_noise_floor_sigma=0)

    def test_readiness_does_not_promote_machine_commit_or_startup(self):
        r=X.readiness()
        for k in ('pending_SpectralMSE_interval_commit_boundary_materialized',
                  'same_tau_drives_exact_real_OU_and_pseudo_cadence',
                  'same_sigma_drives_exact_real_stationary_aw_target',
                  'RS_interval_propagated_through_base_clamp_axis_factors_and_covariance_square',
                  'no_point_RS_representative_selected',
                  'cubic_information_rate_sqrt_not_in_deployed_SpectralMSE_commit'):
            self.assertTrue(r[k])
        for k in ('binary32_commit_arithmetic_correspondence_closed','machine_real_RS_residual_attached_to_commit',
                  'goLive_interval_RS_to_actual_MEKF_commit_closed','source_uniform_complete_startup_reachability_closed',
                  'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
