"""Binary32 shipping sigma-target regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as X
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT


def cfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def sqrt_witness(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)

def exp_witness(x):
    lo,hi,_,_=EXP.enclosure(x); return B.rn32((lo+hi)/2)


class Tests(unittest.TestCase):
    def test_ready_nonstill_sigma_target_uses_noise_subtraction_floor_sqrt_gain_clamp(self):
        c=cfg(); av=B.rn32(F(1,4)); bn=B.rn32(F(1,10)); vn=B.mul(bn,bn)
        pre=max(B.rn32(0),B.sub(av,vn)); vw=max(pre,X.VAR_FLOOR)
        out=X.target(c,var_ready=True,accel_variance=av,band_noise_sigma=bn,
                     still=False,sqrt_result=sqrt_witness(vw))
        self.assertEqual(out.var_noise,vn); self.assertEqual(out.var_wave_pre_attenuation,pre)
        self.assertEqual(out.var_wave,vw); self.assertEqual(out.attenuation,B.rn32(1))
        self.assertEqual(out.scaled_sigma,B.mul(out.sqrt_result,c.sigma_coeff))
        self.assertEqual(out.sigma_target,min(out.scaled_sigma,c.max_sigma))

    def test_still_branch_consumes_same_argument_exp_and_multiplies_before_floor(self):
        c=cfg(); av=B.rn32(F(1,4)); bn=B.rn32(F(1,10)); st=B.rn32(F(1,2))
        vn=B.mul(bn,bn); pre=max(B.rn32(0),B.sub(av,vn)); e=exp_witness(st)
        attenuated=B.mul(pre,e); vw=max(attenuated,X.VAR_FLOOR)
        out=X.target(c,var_ready=True,accel_variance=av,band_noise_sigma=bn,still=True,
                     still_time=st,still_exp_result=e,sqrt_result=sqrt_witness(vw))
        self.assertEqual(out.attenuation,e); self.assertEqual(out.var_wave_attenuated,attenuated)
        self.assertEqual(out.var_wave,vw)

    def test_unready_branch_starts_at_noise_variance_and_applies_sigma_floor(self):
        c=cfg(); bn=B.rn32(F(1,100)); vn=B.mul(bn,bn); vw=X.VAR_FLOOR
        out=X.target(c,var_ready=False,accel_variance=B.rn32(0),band_noise_sigma=bn,
                     still=False,sqrt_result=sqrt_witness(vw))
        self.assertEqual(out.var_total,vn); self.assertEqual(out.var_wave_pre_attenuation,0)
        self.assertGreaterEqual(out.sigma_target,X.SIGMA_FLOOR)

    def test_detached_sqrt_and_still_exp_fail_closed(self):
        c=cfg(); av=B.rn32(F(1,4)); bn=B.rn32(F(1,10)); vn=B.mul(bn,bn)
        pre=max(B.rn32(0),B.sub(av,vn)); vw=max(pre,X.VAR_FLOOR)
        with self.assertRaisesRegex(ValueError,'sqrt witness detached'):
            X.target(c,var_ready=True,accel_variance=av,band_noise_sigma=bn,still=False,
                     sqrt_result=B.rn32(2))
        with self.assertRaisesRegex(ValueError,'stillness exp witness detached'):
            X.target(c,var_ready=True,accel_variance=av,band_noise_sigma=bn,still=True,
                     still_time=B.rn32(F(1,2)),still_exp_result=B.rn32(F(1,2)),
                     sqrt_result=sqrt_witness(vw))

    def test_readiness_keeps_upstream_and_libm_correspondence_open(self):
        r=X.readiness()
        for k in ('shipping_sigma_target_source_shape_matches','noise_variance_subtraction_and_zero_floor_binary32_materialized',
                  'optional_stillness_attenuation_binary32_multiply_materialized','variance_1e_minus6_floor_binary32_materialized',
                  'sigma_sqrt_witness_bound_to_same_rounded_var_wave','sigma_gain_max_clamp_and_unready_floor_binary32_materialized'):
            self.assertTrue(r[k])
        for k in ('stillness_exp_target_libm_correspondence_closed','sigma_sqrt_target_libm_correspondence_closed',
                  'upstream_accel_variance_binary32_production_closed','upstream_band_noise_sigma_binary32_production_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
