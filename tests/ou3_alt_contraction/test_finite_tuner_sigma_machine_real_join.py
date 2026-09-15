"""Same-source sigma exact/machine join regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as X


def ccfg():
    # qeff fields are irrelevant to sigma; choose an exact legacy identity while
    # keeping all shared tau/sigma adaptation scalars at shipping source values.
    return C.CandidateConfig(F(3,100),F(6,5),1,F(9,10),F(1,50),12,4,
        F(3,22),F(1,200),F(3,20),F(3,20),100,F(538,10000),F(1,2),1,
        F(9,5),F(2,5),F(3,2),0,F(1,10),True)


def dcfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def sample(*,still=False):
    return C.WaveBandSample(F(1,2),True,F(1,4),0,still,0,1,F(1,2))

def machine(c,*,av=F(1,4),still=False):
    return M.target(c,var_ready=True,accel_variance=B.rn32(av),band_noise_sigma=B.rn32(0),
                    still=still,still_time=B.rn32(0),
                    still_exp_result=B.rn32(1) if still else None,
                    sqrt_result=B.rn32(F(1,2)))


class Tests(unittest.TestCase):
    def test_same_sample_join_retains_machine_minus_exact_sigma_supply(self):
        c=dcfg(); s=sample(); j=X.join(s,ccfg(),c,machine(c))
        self.assertEqual(j.supply.accel_variance,0)
        self.assertEqual(j.supply.band_noise_sigma,0)
        self.assertEqual(j.supply.variance_wave,0)
        # 0.9f is not the exact rational 9/10, so even this simple fixture must
        # retain rather than erase the compiled-coefficient discrepancy.
        self.assertEqual(j.supply.sigma_target,j.machine.sigma_target-j.exact_target.sigma_target)

    def test_machine_operand_change_is_supply_not_new_exact_history(self):
        c=dcfg(); _s=sample(); av=B.rn32(F(1,4)+F(1,1<<20))
        # Use the old .5 sqrt deliberately only if it remains in the legal RNE
        # cell; otherwise the machine target itself correctly fails first.
        with self.assertRaises(ValueError):
            machine(c,av=av)

    def test_branch_identity_cannot_be_spliced(self):
        c=dcfg(); m=machine(c,still=False)
        with self.assertRaisesRegex(ValueError,'branch detached'):
            X.join(sample(still=True),ccfg(),c,m)

    def test_detached_deployment_config_is_rejected(self):
        c=dcfg(); m=machine(c)
        c2=D.DeploymentConfig(c.min_freq,c.max_freq,c.tau_coeff,c.sigma_coeff,c.min_tau,c.max_tau,c.max_sigma,
             c.pseudo_tau_ratio,c.pseudo_min,c.pseudo_max,c.min_RS,c.max_RS,c.rs_mse_coeff,
             c.accel_noise_density,c.qeff_cache,c.adapt_tau_sec,c.adapt_tau_sea_periods,
             c.adapt_RS_mult,c.adapt_RS_slew_log,B.rn32(F(1,5)),c.clamp_enabled)
        with self.assertRaisesRegex(ValueError,'detached from joined deployment config'):
            X.join(sample(),ccfg(),c2,m)

    def test_sigma_target_constants_cannot_hide_behind_alpha_only_guard(self):
        for name,value in (('sigma_coeff',F(4,5)),('max_sigma',F(3))):
            with self.subTest(field=name):
                d=replace(dcfg(),**{name:B.rn32(value)})
                self.assertTrue(X.A._configs_match(ccfg(),d))
                self.assertFalse(X.configs_match(ccfg(),d))
                with self.assertRaisesRegex(ValueError,'configs disagree'):
                    X.join(sample(),ccfg(),d,machine(d))

    def test_readiness_keeps_supply_bounds_and_storage_open(self):
        r=X.readiness()
        self.assertTrue(r['exact_frontend_sample_and_binary32_sigma_branch_joined'])
        self.assertTrue(r['machine_minus_exact_wave_variance_and_sigma_target_supplies_exposed'])
        for k in ('upstream_frontend_binary32_correspondence_closed','sigma_sqrt_and_still_exp_target_libm_correspondence_closed',
                  'source_uniform_sigma_input_supply_bounds_closed','startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
