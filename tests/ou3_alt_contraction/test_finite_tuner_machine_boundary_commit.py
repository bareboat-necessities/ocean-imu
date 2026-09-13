"""Binary32 coherent machine TuneState boundary regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary_commit as X
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as M
import test_finite_source_bound_live_word as BASE


def cfg(): return BASE.root_state().runtime.commit_cfg


class Tests(unittest.TestCase):
    def test_no_pending_boundary_is_complete_machine_identity(self):
        s=M.initial(); out=X.apply(s,cfg(),live=False)
        self.assertIs(out.state,s); self.assertIsNone(out.separate); self.assertIsNone(out.fma)
        with self.assertRaisesRegex(ValueError,'no-pending'):
            X.apply(s,cfg(),live=False,separate_band_noise_floor_sigma=B.rn32(F(1,20)))

    def test_preLive_pending_commits_tau_period_and_sigma_but_not_RS(self):
        s=replace(M.initial(),pending=True); bn=B.rn32(F(3,50))
        out=X.apply(s,cfg(),live=False,separate_band_noise_floor_sigma=bn,fma_band_noise_floor_sigma=bn)
        self.assertFalse(out.state.pending); self.assertEqual(out.state.tau,s.tau); self.assertEqual(out.state.sigma,s.sigma); self.assertEqual(out.state.rs,s.rs)
        self.assertIsNone(out.separate.rs); self.assertIsNone(out.fma.rs)
        self.assertEqual(out.separate.tau_command,max(X.TAU_FLOOR,s.tau.separate))
        q=X._cfg32(cfg()); requested=B.mul(q['pseudo_tau_ratio'],s.tau.separate)
        self.assertEqual(out.separate.pseudo_requested,requested)
        floor=max(X.SIGMA_FLOOR,bn); self.assertEqual(out.separate.sigma_z,max(floor,s.sigma.separate))

    def test_Live_pending_commits_RS_from_same_mode_snapshot(self):
        s=replace(M.initial(),pending=True); sb=B.rn32(F(1,20)); fb=B.rn32(F(3,50))
        out=X.apply(s,cfg(),live=True,separate_band_noise_floor_sigma=sb,fma_band_noise_floor_sigma=fb)
        self.assertIsNotNone(out.separate.rs); self.assertIsNotNone(out.fma.rs)
        self.assertEqual(out.separate.rs.stored_RS,s.rs.separate); self.assertEqual(out.fma.rs.stored_RS,s.rs.fma)
        self.assertEqual(out.separate.band_noise_floor_sigma,sb); self.assertEqual(out.fma.band_noise_floor_sigma,fb)
        self.assertEqual(out.separate.aw_covariance_diag,tuple(B.mul(v,v) for v in out.separate.aw_std))

    def test_pending_requires_both_mode_specific_band_noise_floor_witnesses(self):
        s=replace(M.initial(),pending=True)
        with self.assertRaisesRegex(ValueError,'each compiler track'):
            X.apply(s,cfg(),live=True,separate_band_noise_floor_sigma=B.rn32(F(1,20)))

    def test_readiness_closes_machine_boundary_not_upstream_or_Eigen(self):
        r=X.readiness()
        for k in ('shipping_pending_common_TuneState_boundary_source_shape_matches','boundary_consumes_no_exact_real_TuneState',
                  'same_mode_tau_sigma_RS_snapshot_committed_together','tau_setter_floor_and_pseudo_period_binary32_graph_materialized',
                  'sigma_floor_horizontal_scale_and_covariance_square_binary32_graph_materialized',
                  'Live_SpectralMSE_RS_commit_composed_from_same_machine_snapshot','preLive_boundary_commits_OU_but_not_RS',
                  'pending_bit_cleared_without_advancing_candidate_ledgers'):
            self.assertTrue(r[k])
        for k in ('upstream_band_noise_floor_binary32_production_closed','Eigen_aw_stationary_std_execution_correspondence_closed',
                  'Eigen_RS_noise_execution_correspondence_closed','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
