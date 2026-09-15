"""SpectralMSE exact-candidate to machine-input ancestry regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SM
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SJ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_input_join as X


def ccfg():
    return C.CandidateConfig(F(3,100),F(6,5),1,F(9,10),F(1,50),12,4,
        F(3,22),F(1,200),F(3,20),F(3,20),100,F(538,10000),F(1,2),1,
        F(9,5),F(2,5),F(3,2),0,F(1,10),True)

def dcfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def exact_sample():
    return C.WaveBandSample(F(1,2),True,F(1,4),0,False,0,1,F(1,2))

def sigma_join():
    d=dcfg(); m=SM.target(d,var_ready=True,accel_variance=B.rn32(F(1,4)),
        band_noise_sigma=B.rn32(0),still=False,still_time=B.rn32(0),
        sqrt_result=B.rn32(F(1,2)))
    return SJ.join(exact_sample(),ccfg(),d,m)

def tau_step():
    _c=ccfg(); f=B.rn32(F(1,2)); target=B.rn32(1); sea=B.rn32(1)
    adapt=B.rn32(F(2,5)); e=B.rn32(1); alpha=B.rn32(0); prev=B.rn32(F(11,10))
    # Structural fixture: the input join consumes the already-qualified TauStep;
    # tau arithmetic qualification itself has independent exhaustive tests.
    return T.TauStep(prev,f,target,sea,adapt,e,alpha,prev,prev)


class Tests(unittest.TestCase):
    def test_machine_target_uses_same_tau_and_sigma_operands(self):
        tj=tau_step(); sj=sigma_join()
        out=X.join(tj,sj,pow_result=B.rn32(1),sqrt_result=B.rn32(F(1,3)))
        self.assertEqual(out.machine_target.tau_target,tj.tau_target)
        self.assertEqual(out.machine_target.sigma_target,sj.machine.sigma_target)
        self.assertEqual(out.machine_target.variance_wave,sj.machine.var_wave)
        self.assertEqual(out.spectral.target,out.machine_target)
        self.assertEqual(out.input_supply.sigma_target,sj.supply.sigma_target)
        self.assertEqual(out.input_supply.variance_wave,sj.supply.variance_wave)

    def test_compiled_sigma_difference_is_not_hidden_in_spectral_libm_residual(self):
        out=X.join(tau_step(),sigma_join(),pow_result=B.rn32(1),sqrt_result=B.rn32(F(1,3)))
        self.assertNotEqual(out.input_supply.sigma_target,0)
        # Local spectral residual is defined at the machine-input target; it is
        # a different proof supply than exact-candidate -> machine-input shift.
        self.assertEqual(out.spectral.exact.target,out.machine_target)
        self.assertNotEqual(out.spectral.exact.target,out.exact_target)

    def test_constructor_rejects_recomputed_sigma_supply(self):
        out=X.join(tau_step(),sigma_join(),pow_result=B.rn32(1),sqrt_result=B.rn32(F(1,3)))
        bad=X.InputSupply(out.input_supply.frequency,out.input_supply.variance_wave,
                          out.input_supply.tau_target,0)
        with self.assertRaisesRegex(ValueError,'sigma input supply detached'):
            X.Join(out.exact_target,out.machine_target,out.tau_step,out.sigma_join,bad,out.spectral)

    def test_readiness_separates_input_and_libm_obligations(self):
        r=X.readiness()
        self.assertTrue(r['machine_input_displacement_kept_separate_from_local_spectral_libm_roundoff'])
        self.assertTrue(r['sigma_variance_and_target_supplies_reused_without_duplication'])
        for k in ('spectral_input_displacement_to_exact_candidate_RS_bound_closed','target_libm_correspondence_closed',
                  'source_uniform_spectral_input_supply_bounds_closed','startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
