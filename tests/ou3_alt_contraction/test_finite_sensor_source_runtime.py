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
    def make(self,D=S.IDENTITY3,acc_res=(0,0,0)):
        state=FC.root('A'); phys=state.reference
        omega=(F(1,50),0,0); ng=(F(1,100000),0,0)
        gint=tuple(omega[i]+phys.gyro_bias[i]+ng[i] for i in range(3))
        inertial=(phys.acceleration[0],phys.acceleration[1],phys.acceleration[2]-G)
        fint=S.q_rotate(phys.q_world_to_body,inertial); aint=tuple(fint[i]+phys.beta[i]+F(acc_res[i]) for i in range(3))
        DT=tuple(zip(*D)); graw=S.mv3(DT,gint); araw=S.mv3(DT,aint)
        return state,S.RawImuSample(phys,omega,ng,acc_res,graw,araw,(0,0,G),D)

    def test_internal_gyro_relation(self):
        state,sample=self.make(); corrected=sample.bias_corrected_internal_gyro(state.z[3:6])
        self.assertEqual(sample.internal_gyro,(F(2011,100000),0,0)); self.assertEqual(corrected,sample.required_bias_corrected_relation(state.z[3:6])); self.assertEqual(corrected,(F(2001,100000),0,0))

    def test_internal_accel_relation(self):
        _,sample=self.make(); self.assertEqual(sample.internal_accel,(F(11,100),0,-G))

    def test_nonzero_deheel_links_distinct_raw_and_internal_packets(self):
        D=((1,0,0),(0,0,1),(0,-1,0)); _,sample=self.make(D)
        self.assertNotEqual(sample.raw_accel_body,sample.internal_accel); self.assertEqual(S.mv3(D,sample.raw_accel_body),sample.internal_accel); self.assertEqual(S.mv3(D,sample.raw_gyro_body),sample.internal_gyro)

    def test_nonorthogonal_deheel_rejected(self):
        with self.assertRaises(ValueError): self.make(((1,0,0),(0,2,0),(0,0,1)))

    def test_temperature_removed_core_observation_and_nu_share_same_raw_residual(self):
        _,sample=self.make(acc_res=(F(3,100),F(-1,100),F(1,50)))
        c=S.AccelConditioning(2,(F(1,100),F(1,50),F(-1,100)))
        out=S.finite_accel_core_observation(sample,c)
        modeled=(F(1,50),F(1,25),F(-1,50))
        self.assertEqual(out.observed,tuple(sample.internal_accel[i]-modeled[i] for i in range(3)))
        self.assertEqual(out.nu_acc,(F(1,100),F(-1,20),F(1,25)))

    def test_finite_acc_bridge_enforces_declared_zero_lever_branch(self):
        _,sample=self.make()
        with self.assertRaises(ValueError): S.finite_accel_core_observation(sample,S.AccelConditioning(0,(0,0,0),(1,0,0)))

    def test_detached_sensor_values_rejected(self):
        state,sample=self.make()
        with self.assertRaises(ValueError): S.RawImuSample(sample.physical,sample.omega_sample_internal,sample.gyro_residual_internal,sample.accel_residual_internal,(0,0,0),sample.raw_accel_body,sample.gravity_world,sample.deheel_body_to_internal)
        bad=list(sample.raw_accel_body); bad[0]+=1
        with self.assertRaises(ValueError): S.RawImuSample(sample.physical,sample.omega_sample_internal,sample.gyro_residual_internal,sample.accel_residual_internal,sample.raw_gyro_body,bad,sample.gravity_world,sample.deheel_body_to_internal)
        with self.assertRaises(ValueError): S.assert_prediction_gyro(sample,state,(0,0,0))
        with self.assertRaises(ValueError): S.assert_acc_measurement_input(sample,(0,0,0))

    def test_private_vertical_consumes_raw_body_pair(self):
        state,_=self.make(); phys=replace(state.reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0))
        sample=S.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-G),(0,0,G))
        out=S.vertical_step_from_raw(V.State(initialized=True),V.Config(0,0,G,20),sample,dt=F(1,200),accel_invnorm=V.InvSqrtWitness(G*G,1/G),quat_invnorm=V.InvSqrtWitness(1,1)); self.assertEqual(out.vertical_accel,0)

    def test_readiness_fail_closed(self):
        r=S.readiness(); self.assertTrue(r['raw_body_to_internal_deheel_map_materialized']); self.assertTrue(r['raw_accel_to_temperature_removed_core_observation_bridge']); self.assertTrue(r['nu_acc_from_same_raw_residual_and_temperature_model']); self.assertTrue(r['zero_lever_theorem_branch_enforced']); self.assertFalse(r['temperature_and_k_a_runtime_ancestry_attached']); self.assertFalse(r['sensor_residual_source_bounds_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
