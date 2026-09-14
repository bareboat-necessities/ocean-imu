"""Binary32 SpectralMSE R_S commit regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_rs_commit_binary32 as X


class Tests(unittest.TestCase):
    def test_default_spectral_commit_clamps_scales_and_squares(self):
        r=B.rn32(F(13,10)); lo=B.rn32(F(15,100)); hi=B.rn32(100)
        out=X.commit(r,min_RS=lo,max_RS=hi)
        self.assertEqual(out.RSbase,r)
        self.assertEqual(out.RSb,r)
        self.assertEqual(out.rs_z,r)
        self.assertEqual(out.sigma_x,B.mul(r,X.DEFAULT_X_FACTOR))
        self.assertEqual(out.sigma_y,B.mul(r,X.DEFAULT_Y_FACTOR))
        self.assertEqual(out.sigma_z,r)
        self.assertEqual(out.covariance_diag,(B.mul(out.sigma_x,out.sigma_x),
                                              B.mul(out.sigma_y,out.sigma_y),
                                              B.mul(out.sigma_z,out.sigma_z)))

    def test_min_and_max_clamps_are_literal(self):
        lo=B.rn32(F(15,100)); hi=B.rn32(100)
        a=X.commit(B.rn32(F(1,100)),min_RS=lo,max_RS=hi)
        b=X.commit(B.rn32(128),min_RS=lo,max_RS=hi)
        self.assertEqual(a.RSbase,lo); self.assertEqual(b.RSbase,hi)

    def test_source_qualified_rs_scale_is_applied_before_axis_factors(self):
        lo=B.rn32(F(15,100)); hi=B.rn32(100); scale=B.rn32(F(3,4)); r=B.rn32(2)
        out=X.commit(r,min_RS=lo,max_RS=hi,rs_scale=scale)
        self.assertEqual(out.rs_z,B.mul(r,scale))
        self.assertEqual(out.sigma_x,B.mul(out.rs_z,X.DEFAULT_X_FACTOR))

    def test_invalid_scale_does_not_sneak_through_strong_entry(self):
        with self.assertRaisesRegex(ValueError,'source-qualified'):
            X.commit(B.rn32(1),min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100),rs_scale=B.rn32(0))

    def test_readiness_keeps_eigen_and_master_word_open(self):
        r=X.readiness()
        for k in ('shipping_SpectralMSE_RS_commit_source_shape_matches',
                  'SpectralMSE_information_rate_scale_is_literal_one',
                  'stored_RS_clamp_binary32_materialized',
                  'rs_scale_and_anisotropic_factor_binary32_multiplies_materialized',
                  'per_axis_covariance_square_binary32_graph_materialized'):
            self.assertTrue(r[k])
        for k in ('Eigen_set_RS_noise_execution_correspondence_closed',
                  'compiler_vectorization_rounding_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
