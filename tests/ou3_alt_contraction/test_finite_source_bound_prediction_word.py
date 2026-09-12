"""Prediction-root provenance regressions; not storage/contraction evidence."""
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as QX
import test_finite_source_bound_live_word as BASE

PASS=QX.PSDWitness(True)
ROOT=Path(__file__).resolve().parents[2]
CORE=ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
WRAPPER=ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'


def operands(state):
    witness,segment,raw,dynamic=BASE.next_imu_operands(state)
    for key in ('angular','ou','qaxis','bias'):
        dynamic.pop(key,None)
    return witness,segment,raw,dynamic


def root_args():
    return dict(ou_alpha=F(199,200),ou_em1=F(-1,200),
                qaxis_marginal_psd=(PASS,PASS,PASS),
                qaxis_final_psd=(PASS,PASS,PASS))


class Tests(unittest.TestCase):
    def test_prediction_roots_are_derived_from_same_source_and_active_state(self):
        s=BASE.root_state(); witness,segment,raw,_=operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        r=X.build(s,physical,raw,**root_args())
        core=s.live.live.live.mekf; active=s.live.live.live.active
        self.assertEqual(r.angular.w,raw.required_bias_corrected_relation(core.z[3:6]))
        self.assertEqual(r.angular.h,segment.h)
        self.assertEqual(r.ou.h,segment.h)
        self.assertEqual(r.ou.tau,active.tau)
        self.assertEqual(r.ou.alpha,F(199,200)); self.assertEqual(r.ou.em1,F(-1,200))
        self.assertEqual(r.qaxis.sigma_aw,active.Sigma_aw)
        self.assertFalse(r.qaxis.correlated)
        self.assertEqual(r.qaxis.machine_epsilon,F(1,1 << 23))
        self.assertFalse(r.bias.active)
        self.assertEqual(r.bias.tau_b,F(5000))
        self.assertEqual(r.bias.phi_b,F(1)); self.assertEqual(r.bias.em1_2,F(0))
        self.assertEqual(r.bias.Q_bacc,X.SHIPPING_BA_Q)

    def test_source_bound_entry_executes_without_free_prediction_roots(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        out=X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='imu-1',
                       **root_args(),**dynamic)
        self.assertEqual(out.state.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live.live.live.mekf.reference,segment.after)

    def test_wrong_source_or_ordinal_fails_before_runtime_execution(self):
        s=BASE.root_state(); _,segment,raw,dynamic=operands(s)
        bad=SOURCE.StepWitness(2,'root','c2','p0','p2')
        with self.assertRaisesRegex(ValueError,'next source ordinal'):
            X.imu_step(s,witness=bad,segment=segment,raw=raw,packet_id='bad',
                       **root_args(),**dynamic)

    def test_prediction_roots_cannot_be_overridden_by_theorem_caller(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        for key,value in (('ou','detached'),('bias','detached'),('machine_epsilon',F(1,10**7))):
            bad=dict(dynamic); bad[key]=value
            with self.assertRaisesRegex(TypeError,'cannot be overridden'):
                X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='bad',
                           **root_args(),**bad)

    def test_H18_consumes_no_BA_transcendentals_and_A21_requires_both(self):
        held=X._bias_root(SimpleNamespace(mode='H'))
        self.assertFalse(held.active); self.assertEqual(held.phi_b,1); self.assertEqual(held.em1_2,0)
        with self.assertRaisesRegex(ValueError,'consumes no exp/expm1'):
            X._bias_root(SimpleNamespace(mode='H'),bias_phi=F(9,10),bias_em1_2=F(-1,5))
        with self.assertRaisesRegex(ValueError,'requires exp and expm1'):
            X._bias_root(SimpleNamespace(mode='A'),bias_phi=F(999999,1000000))
        active=X._bias_root(SimpleNamespace(mode='A'),bias_phi=F(999999,1000000),
                            bias_em1_2=F(-1,500000))
        self.assertTrue(active.active)
        self.assertEqual(active.tau_b,F(5000))
        self.assertEqual(active.phi_b,F(999999,1000000))
        self.assertEqual(active.em1_2,F(-1,500000))
        self.assertEqual(active.Q_bacc,X.SHIPPING_BA_Q)

    def test_shipping_BA_defaults_and_wrapper_ownership_are_source_locked(self):
        self.assertEqual(X.SHIPPING_BA_TAU,F(5000))
        self.assertEqual(X.SHIPPING_BA_Q,
            ((F(1,4_000_000),0,0),(0,F(1,4_000_000),0),(0,0,F(1,4_000_000))))
        core=CORE.read_text(); wrapper=WRAPPER.read_text()
        self.assertEqual(core.count('Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);'),1)
        self.assertEqual(core.count('T tau_bacc_ = T(5000.0);'),1)
        self.assertIn('const T tau_b = std::max(T(1e-3), tau_bacc_);',core)
        self.assertIn('const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);',core)
        self.assertIn('const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);',core)
        for setter in ('set_acc_bias_time_constant(', 'set_Q_bacc_rw(',
                       'set_acc_bias_ou_stationary_std('):
            self.assertNotIn(setter,wrapper)

    def test_shipping_float_instantiation_locks_qaxis_machine_epsilon(self):
        wrapper=WRAPPER.read_text()
        self.assertGreaterEqual(wrapper.count('Kalman3D_Wave_OU_III<float>'),3)
        self.assertEqual(X.SHIPPING_FLOAT_EPSILON,F(1,1 << 23))

    def test_readiness_advances_only_prediction_root_provenance(self):
        r=X.readiness()
        self.assertTrue(r['prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation'])
        self.assertTrue(r['OU_exp_and_expm1_shipping_results_retained_separately'])
        self.assertTrue(r['Qaxis_machine_epsilon_bound_to_shipping_binary32'])
        self.assertTrue(r['shipping_BA_tau_Q_defaults_bound_at_prediction_entry'])
        self.assertTrue(r['shipping_BA_hold_active_branch_derived_from_current_MEKF_mode'])
        self.assertTrue(r['held_H18_BA_phi_exactly_one_without_transcendental_witness'])
        self.assertTrue(r['active_A21_BA_exp_and_expm1_results_retained_separately'])
        self.assertTrue(r['accelerometer_bias_prediction_root_bound_here'])
        self.assertTrue(r['caller_cannot_override_prediction_roots'])
        self.assertFalse(r['OU_exp_expm1_binary32_relation_closed'])
        self.assertFalse(r['BA_exp_expm1_binary32_relation_closed'])
        self.assertFalse(r['Qaxis_PSD_regularization_deployment_closed'])
        self.assertFalse(r['all_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
