from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as S
from tools.stability.ou3_alt_contraction import finite_startup_disturbance_obstruction as OBSTRUCTION
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_frontend_uniform_bounds as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
import test_finite_startup_joined_machine_history as BASE


# The inherited rational-root fixture is a component, not the shipping source.
# A bounded (commissioned) profile attaches the configured candidate, commit and
# frontend uniform supplies, so that regime has to be entered from the actual
# compiled deployment constants rather than the component scalars.
SHIPPING_TUNER_FIELDS=('min_freq','max_freq','tau_coeff','min_tau','max_tau',
                       'adapt_tau_sec','adapt_tau_sea_periods','sigma_coeff','max_sigma')


def shipping_deployment(cfg):
    return replace(cfg,sigma_coeff=D.SIGMA_COEFF,max_sigma=D.MAX_SIGMA)


def shipping_runtime(runtime,deployment_cfg,**overrides):
    candidate=replace(runtime.candidate_cfg,
        **{n:F(getattr(deployment_cfg,n)) for n in SHIPPING_TUNER_FIELDS},
        clamp_enabled=deployment_cfg.clamp_enabled)
    commit=replace(runtime.commit_cfg,tau_scaled_cadence=True,cubic_rs_law=False,
        pseudo_tau_ratio=D.PSEUDO_RATIO,pseudo_period_min=D.PSEUDO_MIN,
        pseudo_period_max=D.PSEUDO_MAX,pseudo_fixed_period=D.PSEUDO_NOMINAL,
        min_R_S=D.MIN_RS,max_R_S=D.MAX_RS,S_factor=F(1),
        R_S_x_factor=B.rn32(F(72,100)),R_S_y_factor=B.rn32(F(72,100)))
    return replace(runtime,candidate_cfg=candidate,commit_cfg=commit,
        bench_noise_sigma=FRONT.BENCH_SIGMA,boundary_bench_noise_sigma=FRONT.BENCH_SIGMA,
        **overrides)


def history(raw,name='BMI270_COMMISSIONED_V1'):
    r=raw.physical
    return S.History(S.profile(name),r.history_id,r.bias_root,r.bias_family,'gyro-history','acc-history')


class Tests(unittest.TestCase):
    def test_all_bias_families_and_adversarial_packet_bounds(self):
        for name in S.domain()['profiles']:
            p=S.profile(name)
            for family in ('BIAS0','BIAS1','BIAS2'):
                r=replace(OBSTRUCTION.quiet_packet(bias_family=family).physical,
                          acceleration=(0,0,S.A),beta=(0,0,S.BIAS_COMPONENT))
                residual=(0,0,p.accel_residual)
                acc=(0,0,S.A-S.G+S.BIAS_COMPONENT+p.accel_residual)
                raw=S.SENSOR.RawImuSample(r,(0,0,0),(0,0,0),residual,(0,0,0),acc)
                cert=S.check_packet(history(raw,name),raw,ordinal=1)
                seed=SEED.seed(tuple(S.B.rn32(x) for x in acc))
                self.assertLessEqual(cert['computed_norm_lower'],seed.norm)
                self.assertLessEqual(seed.norm,cert['computed_norm_upper'])
                self.assertNotEqual(seed.branch,'acc-norm-too-small')

    def test_constant_residual_obstruction_excluded_by_selected_caps(self):
        raw=OBSTRUCTION.quiet_packet()
        for name in S.domain()['profiles']:
            with self.assertRaisesRegex(ValueError,'accelerometer residual exceeds'):
                S.check_packet(history(raw,name),raw,ordinal=1)

    def test_defined_mahony_projection_supplies_WPE_without_tilt_accuracy(self):
        for name in S.domain()['profiles']:
            r=S.vertical_supply_certificate(S.profile(name))
            self.assertTrue(r['all_nonnegative_finite_norm_words_including_zero_subnormal_covered'])
            self.assertLess(r['computed_quaternion_norm2_upper'],F(139,125))
            self.assertLess(r['vertical_abs_upper_after_defined_Mahony_update'],32)
            self.assertFalse(r['preceding_Mahony_operations_source_uniformly_total'])

    def test_joint_ancestry_clock_and_profile_cannot_be_replaced(self):
        base=BASE.initial(); raw,_=BASE.cold_operands(base); h=history(raw)
        S.check_packet(h,raw,ordinal=1)
        with self.assertRaisesRegex(ValueError,'5ms source clock'):
            S.check_packet(h,raw,ordinal=2)
        with self.assertRaisesRegex(ValueError,'physical/BIAS history'):
            S.check_packet(replace(h,bias_root='different'),raw,ordinal=1)
        with self.assertRaisesRegex(ValueError,'declared commissioned sensor bounds'):
            replace(h.profile,accel_residual=F(10))
        with self.assertRaisesRegex(ValueError,'observer gravity detached'):
            BASE.X.initial(replace(base.runtime,vertical_cfg=replace(base.runtime.vertical_cfg,gravity=F(100))),
                           base.deployment_cfg,sensor_history=h)

    def test_two_actual_joined_steps_and_pending_boundary_keep_profile(self):
        old=BASE.initial(); raw,_=BASE.cold_operands(old); h=history(raw)
        s=BASE.X.initial(old.runtime,old.deployment_cfg,sensor_history=h)
        for _ in range(2):
            raw,kw=BASE.cold_operands(s)
            out=BASE.X.step(s,raw,**kw)
            self.assertIs(out.state.sensor_history,h)
            self.assertTrue(out.separate_source.mahony.vertical.state.initialized)
            self.assertEqual(out.guard.conditioned_acc,kw['machine_acc_body'])
            s=out.state
        held,_=BASE.X.boundary(s)
        self.assertIs(held.sensor_history,h)

    def test_source_implications_do_not_promote_hardware_or_temporal_membership(self):
        r=S.build(); self.assertEqual(S.validate(r),[])
        self.assertTrue(r['source_uniform_seed_norm_margin_closed'])
        for key in ('temporal_budget_follows_from_peak_checks','old_Mahony_invariant_covers_new_profiles',
                    'hardware_admission_qualified','target_compiler_libm_qualified','storage_search_allowed'):
            self.assertFalse(r[key])
        r['old_Mahony_invariant_covers_new_profiles']=True
        self.assertTrue(S.validate(r))

    def test_full_startup_wpe_product_retains_selected_profile(self):
        from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as X
        from test_finite_startup_wpe_machine_history import mode_witness
        from test_finite_wpe_uniform_bounds import config
        old=BASE.initial(); raw,_=BASE.cold_operands(old); h=history(raw)
        deployment=shipping_deployment(old.deployment_cfg)
        runtime=shipping_runtime(old.runtime,deployment,wpe_cfg=config())
        s=X.initial(runtime,deployment,sensor_history=h)
        for _ in range(2):
            raw,kw=BASE.cold_operands(s.base)
            guard=BASE.X.GUARD.step(s.base.guard,s.base.guard_cfg,dt=kw['machine_dt'],
                raw_gyro=kw['machine_gyro_body'],raw_acc=kw['machine_acc_body'],**kw['guard_witnesses'])
            src,_=BASE.VSBASE.source_witness(s.base.separate_source,guard,s.base.runtime,kw['dt'])
            w=mode_witness(s.wpe,src.band_input,kw['dt'])
            s=X.step(s,raw,separate_wpe=w,fma_wpe=w,**kw).state
        held,_=X.boundary(s)
        self.assertIs(held.base.sensor_history,h)
        self.assertIs(held.wpe,s.wpe)
        with self.assertRaisesRegex(ValueError,'bounds profile detached'):
            replace(s,base=replace(s.base,sensor_history=None))


if __name__=='__main__': unittest.main()
