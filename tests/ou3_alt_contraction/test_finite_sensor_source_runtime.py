"""Raw sensor/source provenance regressions; no source-bound promotion."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as S
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
import test_finite_core as FC

G=F(196133,20000)


class Tests(unittest.TestCase):
    def make(self):
        state=FC.root('A')
        phys=state.reference
        omega=(F(1,50),0,0); ng=(F(1,100000),0,0)
        gyro=tuple(omega[i]+phys.gyro_bias[i]+ng[i] for i in range(3))
        inertial=(phys.acceleration[0],phys.acceleration[1],phys.acceleration[2]-G)
        # identity physical quaternion in FC.root, raw residual zero.
        acc=tuple(inertial[i]+phys.beta[i] for i in range(3))
        return state,S.RawImuSample(phys,omega,ng,(0,0,0),gyro,acc,(0,0,G))

    def test_raw_gyro_is_physical_rate_plus_true_bias_plus_same_residual(self):
        state,sample=self.make()
        self.assertEqual(sample.raw_gyro,(F(2011,100000),0,0))
        corrected=sample.bias_corrected_gyro(state.z[3:6])
        self.assertEqual(corrected,sample.required_bias_corrected_relation(state.z[3:6]))
        self.assertEqual(corrected,(F(2001,100000),0,0))

    def test_raw_accel_is_rotated_specific_force_plus_true_beta_plus_residual(self):
        _,sample=self.make()
        self.assertEqual(sample.raw_accel,(F(11,100),0,-G))

    def test_detached_raw_sensor_values_are_rejected(self):
        state,sample=self.make()
        with self.assertRaises(ValueError):
            S.RawImuSample(sample.physical,sample.omega_sample_body,sample.gyro_residual,
                           sample.accel_residual_raw,(0,0,0),sample.raw_accel,sample.gravity_world)
        bad=list(sample.raw_accel); bad[0]+=1
        with self.assertRaises(ValueError):
            S.RawImuSample(sample.physical,sample.omega_sample_body,sample.gyro_residual,
                           sample.accel_residual_raw,sample.raw_gyro,bad,sample.gravity_world)
        with self.assertRaises(ValueError): S.assert_prediction_gyro(sample,state,(0,0,0))
        with self.assertRaises(ValueError): S.assert_acc_measurement_input(sample,(0,0,0))

    def test_private_vertical_consumes_exact_same_raw_pair(self):
        state,_=self.make()
        phys=replace(state.reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0))
        sample=S.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-G),(0,0,G))
        out=S.vertical_step_from_raw(V.State(initialized=True),V.Config(0,0,G,20),sample,
             dt=F(1,200),accel_invnorm=V.InvSqrtWitness(G*G,1/G),quat_invnorm=V.InvSqrtWitness(1,1))
        self.assertEqual(out.vertical_accel,0)
        with self.assertRaises(TypeError):
            S.vertical_step_from_raw(V.State(initialized=True),V.Config(0,0,G,20),sample,
              dt=F(1,200),gyro=(0,0,0),accel_invnorm=V.InvSqrtWitness(G*G,1/G),quat_invnorm=V.InvSqrtWitness(1,1))

    def test_projective_quaternion_rotation_is_scale_invariant(self):
        v=(1,2,3)
        self.assertEqual(S.q_rotate((1,0,0,0),v),S.q_rotate((2,0,0,0),v))

    def test_readiness_keeps_sensor_bounds_and_conditioning_open(self):
        r=S.readiness()
        self.assertTrue(r['raw_gyro_physical_bias_residual_identity'])
        self.assertTrue(r['omega_hat_equals_omega_sample_plus_e_bg_plus_n_g'])
        self.assertFalse(r['raw_accel_residual_to_nu_acc_temperature_conditioning_bridge'])
        self.assertFalse(r['sensor_residual_source_bounds_attached'])
        self.assertFalse(r['binary32_sensor_conversion_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
