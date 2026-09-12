"""Prediction/model-root provenance regressions; not storage evidence."""
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as QX
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
import test_finite_source_bound_live_word as BASE

PASS=QX.PSDWitness(True)
ROOT=Path(__file__).resolve().parents[2]
CORE=ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
WRAPPER=ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'


def operands(state):
    witness,segment,raw,dynamic=BASE.next_imu_operands(state)
    for key in ('angular','ou','qaxis','bias','accel_conditioning'):
        dynamic.pop(key,None)
    return witness,segment,raw,dynamic


def root_args():
    # For the fixture tau=1 and h=1/200.  These values are the rigorous
    # alternating-series lower bounds for exp(-x) and expm1(-x), not claims
    # about the target libm's rounded outputs.
    return dict(temperature_c=F(35),ou_alpha=F(199,200),ou_em1=F(-1,200),
                qaxis_marginal_psd=(PASS,PASS,PASS),
                qaxis_final_psd=(PASS,PASS,PASS))


class Tests(unittest.TestCase):
    def test_prediction_roots_are_derived_from_same_source_and_active_state(self):
        s=BASE.root_state(); witness,segment,raw,_=operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        args=root_args(); args.pop('temperature_c')
        r=X.build(s,physical,raw,**args)
        core=s.live.live.live.mekf; active=s.live.live.live.active
        self.assertEqual(r.angular.w,raw.required_bias_corrected_relation(core.z[3:6]))
        self.assertEqual(r.angular.h,segment.h)
        self.assertEqual(r.ou.h,segment.h); self.assertEqual(r.ou.tau,active.tau)
        self.assertEqual(r.ou.alpha,F(199,200)); self.assertEqual(r.ou.em1,F(-1,200))
        self.assertEqual(r.qaxis.sigma_aw,active.Sigma_aw)
        self.assertFalse(r.qaxis.correlated)
        self.assertEqual(r.qaxis.machine_epsilon,F(1,1 << 23))
        self.assertFalse(r.bias.active)
        self.assertEqual(r.bias.tau_b,F(5000)); self.assertEqual(r.bias.phi_b,F(1))
        self.assertEqual(r.bias.em1_2,F(0)); self.assertEqual(r.bias.Q_bacc,X.SHIPPING_BA_Q)

    def test_source_bound_entry_executes_without_free_prediction_or_accel_model_roots(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        out=X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='imu-1',
                       **root_args(),**dynamic)
        self.assertEqual(out.state.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live.live.live.mekf.reference,segment.after)

    def test_detached_OU_transcendental_roots_are_rejected(self):
        s=BASE.root_state(); witness,segment,raw,_=operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        args=root_args(); args.pop('temperature_c')
        bad=dict(args); bad['ou_alpha']=F(9,10)
        with self.assertRaisesRegex(ValueError,'exp root detached'):
            X.build(s,physical,raw,**bad)
        bad=dict(args); bad['ou_em1']=F(-1,10)
        with self.assertRaisesRegex(ValueError,'expm1 root detached'):
            X.build(s,physical,raw,**bad)

    def test_temperature_input_builds_shipping_conditioning_not_free_coefficient(self):
        c=X._accel_conditioning(F(37))
        self.assertIsInstance(c,SENSOR.AccelConditioning)
        self.assertEqual(c.temperature_delta,F(2))
        self.assertEqual(c.k_a_hat_internal,X.SHIPPING_KA)
        self.assertEqual(c.lever_internal,(0,0,0))

    def test_wrong_source_or_ordinal_fails_before_runtime_execution(self):
        s=BASE.root_state(); _,segment,raw,dynamic=operands(s)
        bad=SOURCE.StepWitness(2,'root','c2','p0','p2')
        with self.assertRaisesRegex(ValueError,'next source ordinal'):
            X.imu_step(s,witness=bad,segment=segment,raw=raw,packet_id='bad',
                       **root_args(),**dynamic)

    def test_model_roots_cannot_be_overridden_by_theorem_caller(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        detached=SENSOR.AccelConditioning(0,(0,0,0))
        for key,value in (('ou','detached'),('bias','detached'),
                          ('machine_epsilon',F(1,10**7)),('accel_conditioning',detached)):
            bad=dict(dynamic); bad[key]=value
            with self.assertRaisesRegex(TypeError,'cannot be overridden'):
                X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='bad',
                           **root_args(),**bad)

    def test_H18_consumes_no_BA_transcendentals_and_A21_requires_both(self):
        h=F(1,200)
        held=X._bias_root(SimpleNamespace(mode='H'),h=h)
        self.assertFalse(held.active); self.assertEqual(held.phi_b,1); self.assertEqual(held.em1_2,0)
        with self.assertRaisesRegex(ValueError,'consumes no exp/expm1'):
            X._bias_root(SimpleNamespace(mode='H'),h=h,bias_phi=F(9,10),bias_em1_2=F(-1,5))
        with self.assertRaisesRegex(ValueError,'requires h plus exp and expm1'):
            X._bias_root(SimpleNamespace(mode='A'),h=h,bias_phi=F(999999,1000000))
        # The fixture uses conservative real enclosure endpoints around the
        # distinct h/tau_b and 2h/tau_b shipping calls.
        active=X._bias_root(SimpleNamespace(mode='A'),h=h,
                            bias_phi=F(999999,1000000),bias_em1_2=F(-1,500000))
        self.assertTrue(active.active); self.assertEqual(active.tau_b,F(5000))
        self.assertEqual(active.phi_b,F(999999,1000000)); self.assertEqual(active.em1_2,F(-1,500000))
        self.assertEqual(active.Q_bacc,X.SHIPPING_BA_Q)
        with self.assertRaisesRegex(ValueError,'BA exp root detached'):
            X._bias_root(SimpleNamespace(mode='A'),h=h,bias_phi=F(99,100),bias_em1_2=F(-1,500000))

    def test_shipping_BA_temperature_and_zero_lever_defaults_are_source_locked(self):
        self.assertEqual(X.SHIPPING_BA_TAU,F(5000))
        self.assertEqual(X.SHIPPING_BA_Q,
            ((F(1,4_000_000),0,0),(0,F(1,4_000_000),0),(0,0,F(1,4_000_000))))
        self.assertEqual(X.SHIPPING_TEMP_REF_C,F(35)); self.assertEqual(X.SHIPPING_KA,(F(1,500),)*3)
        core=CORE.read_text(); wrapper=WRAPPER.read_text()
        self.assertEqual(core.count('Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);'),1)
        self.assertEqual(core.count('T tau_bacc_ = T(5000.0);'),1)
        self.assertEqual(core.count('static constexpr T tempC_ref = T(35.0);'),1)
        self.assertEqual(core.count('Vector3 k_a_ = Vector3::Constant(T(0.002));'),1)
        self.assertEqual(core.count('bool   use_imu_lever_arm_       = false;'),1)
        self.assertEqual(core.count('Vector3 r_imu_wrt_cog_body_phys_ = Vector3::Zero();'),1)
        self.assertIn('const T tau_b = std::max(T(1e-3), tau_bacc_);',core)
        self.assertIn('const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);',core)
        self.assertIn('const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);',core)
        for setter in ('set_acc_bias_time_constant(', 'set_Q_bacc_rw(',
                       'set_acc_bias_ou_stationary_std(','set_accel_bias_temp_coeff(',
                       'set_imu_lever_arm_body('):
            self.assertNotIn(setter,wrapper)

    def test_shipping_float_instantiation_locks_qaxis_machine_epsilon(self):
        wrapper=WRAPPER.read_text()
        self.assertGreaterEqual(wrapper.count('Kalman3D_Wave_OU_III<float>'),3)
        self.assertEqual(X.SHIPPING_FLOAT_EPSILON,F(1,1 << 23))

    def test_readiness_advances_only_model_root_provenance(self):
        r=X.readiness()
        for key in ('prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation',
                    'OU_exp_and_expm1_shipping_results_retained_separately',
                    'OU_exp_expm1_real_enclosed_at_same_source_argument',
                    'Qaxis_machine_epsilon_bound_to_shipping_binary32',
                    'shipping_BA_tau_Q_defaults_bound_at_prediction_entry',
                    'shipping_BA_hold_active_branch_derived_from_current_MEKF_mode',
                    'held_H18_BA_phi_exactly_one_without_transcendental_witness',
                    'active_A21_BA_exp_and_expm1_results_retained_separately',
                    'active_A21_BA_exp_expm1_real_enclosed_at_literal_arguments',
                    'accelerometer_temperature_coefficient_bound_to_shipping_default',
                    'zero_lever_accelerometer_scope_enforced_at_source_bound_entry',
                    'temperature_is_explicit_per_sample_input_not_free_model_coefficient',
                    'accelerometer_bias_prediction_root_bound_here',
                    'caller_cannot_override_prediction_or_accel_model_roots'):
            self.assertTrue(r[key])
        for key in ('temperature_history_admissibility_attached',
                    'OU_exp_expm1_binary32_relation_closed','BA_exp_expm1_binary32_relation_closed',
                    'Qaxis_PSD_regularization_deployment_closed',
                    'all_estimator_coefficients_bound_to_same_source_continuation',
                    'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[key])


if __name__=='__main__': unittest.main()
