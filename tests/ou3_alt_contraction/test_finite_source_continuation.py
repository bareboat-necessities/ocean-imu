from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import finite_source_continuation as X
from tools.stability.ou3_alt_contraction import finite_physical_prediction as P
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as S

BIAS0=next(c for c in BIAS.contracts() if c.name=='BIAS0')
VALID_PHI=F.from_float(BIAS0.phi_true.lo)
G=F(980665,100000)


def kin(t, *, live=F(0), beta=(0,0,0)):
    return P.PhysicalKinematics(F(t),(1,0,0,0),(0,0,0),(0,0,0),(0,0,0),
                                (0,0,0),(0,0,0),beta,F(live))


def seg(k, *, phi=VALID_PHI, driver=(0,0,0), beta0=(0,0,0), beta1=None):
    t0=F(k-1,200); t1=F(k,200)
    if beta1 is None:
        beta1=tuple(F(phi)*F(beta0[i])+F(driver[i]) for i in range(3))
    return P.PhysicalSegment(kin(t0,beta=beta0),kin(t1,beta=beta1),
                             (0,0,0),(0,0,0),(0,0,0),F(phi),driver)


def root():
    return X.SourceRoot('hist','generator',0,'BIAS0',
                        'BIAS0:one-physical-history:phi-root-and-driver')


def sensor_root(r): return X.SensorDisturbanceRoot(r,'gyro-noise-history','accel-noise-history')


def wit(k,parent,child,pin,pout):
    return X.StepWitness(k,parent,child,pin,pout)


def qseg(r,k=1):
    parent='root' if k==1 else f'c{k-1}'
    return X.QualifiedPhysicalSegment(r,wit(k,parent,f'c{k}',f'p{k-1}',f'p{k}'),seg(k))


def raw_for(q, *, gyro_res=(0,0,0), accel_res=(0,0,0)):
    p=q.segment.before
    gyro=tuple(p.gyro_bias[i]+F(gyro_res[i]) for i in range(3))
    acc=(F(accel_res[0]),F(accel_res[1]),-G+F(accel_res[2]))
    return S.RawImuSample(p,(0,0,0),gyro_res,accel_res,gyro,acc)


class Tests(unittest.TestCase):
    def test_consecutive_segments_preserve_one_source_and_bias_root(self):
        r=root(); c=X.begin(r)
        c=X.append(c,witness=wit(1,'root','c1','p0','p1'),segment=seg(1))
        c=X.append(c,witness=wit(2,'c1','c2','p1','p2'),segment=seg(2))
        self.assertEqual(c.next_ordinal,3); self.assertFalse(c.complete)
        self.assertEqual(c.steps[0].root,c.steps[1].root)
        self.assertEqual(c.steps[0].segment.after,c.steps[1].segment.before)

    def test_bias_family_phi_is_not_free_per_segment(self):
        r=root(); self.assertLess(BIAS0.phi_true.hi,1.0)
        with self.assertRaises(ValueError):
            X.QualifiedPhysicalSegment(r,wit(1,'root','c1','p0','p1'),seg(1,phi=F(1)))

    def test_source_cell_and_primitive_chains_are_hard(self):
        r=root(); c=X.begin(r)
        c=X.append(c,witness=wit(1,'root','c1','p0','p1'),segment=seg(1))
        with self.assertRaises(ValueError):
            X.append(c,witness=wit(2,'wrong','c2','p1','p2'),segment=seg(2))
        with self.assertRaises(ValueError):
            X.append(c,witness=wit(2,'c1','c2','wrong','p2'),segment=seg(2))

    def test_physical_endpoint_chain_cannot_jump(self):
        r=root(); q1=X.QualifiedPhysicalSegment(r,wit(1,'root','c1','p0','p1'),seg(1))
        b=P.PhysicalKinematics(F(1,200),(0,1,0,0),(0,0,0),(0,0,0),(0,0,0),
                               (0,0,0),(0,0,0),(0,0,0),0)
        a=P.PhysicalKinematics(F(2,200),(0,1,0,0),(0,0,0),(0,0,0),(0,0,0),
                               (0,0,0),(0,0,0),(0,0,0),0)
        s2=P.PhysicalSegment(b,a,(0,0,0),(0,0,0),(0,0,0),VALID_PHI,(0,0,0))
        q2=X.QualifiedPhysicalSegment(r,wit(2,'c1','c2','p1','p2'),s2)
        with self.assertRaises(ValueError): X.Continuation(r,(q1,q2))

    def test_raw_packet_is_owned_by_same_physical_and_disturbance_root(self):
        r=root(); q=qseg(r); sr=sensor_root(r)
        raw=raw_for(q,gyro_res=(F(1,10000),0,0),accel_res=(F(1,10),0,0))
        packet=X.qualify_raw_imu(q,sr,raw,'imu-1')
        self.assertEqual(packet.raw,raw); self.assertEqual(packet.physical,q)
        self.assertEqual(packet.sensor_root.gyro_residual_history_id,'gyro-noise-history')
        self.assertEqual(packet.sensor_root.accel_residual_history_id,'accel-noise-history')

    def test_raw_packet_or_sensor_root_cannot_restart(self):
        r=root(); q1=qseg(r,1); q2=qseg(r,2); raw2=raw_for(q2)
        sr=sensor_root(r)
        with self.assertRaises(ValueError): X.qualify_raw_imu(q1,sr,raw2,'imu-wrong-endpoint')
        other=X.SourceRoot('other','generator',0,'BIAS0',
                           'BIAS0:one-physical-history:phi-root-and-driver')
        with self.assertRaises(ValueError):
            X.qualify_raw_imu(q1,X.SensorDisturbanceRoot(other,'g','a'),raw_for(q1),'imu-wrong-root')


    def test_fresh_origin_is_checked_without_fabricating_transition(self):
        r=root(); p=kin(0)
        origin=X.origin_endpoint(r,p)
        self.assertEqual(origin.endpoint,p); self.assertEqual(origin.root,r)
        with self.assertRaisesRegex(ValueError,'one-time Live origin'):
            X.origin_endpoint(r,kin(F(1,200)))
        with self.assertRaisesRegex(ValueError,'acceleration vector cap'):
            X.origin_endpoint(r,replace(p,acceleration=(100,0,0)))

    def test_async_endpoint_is_not_a_free_matching_state(self):
        r=root(); q=qseg(r)
        before=X.endpoint(q,'before'); after=X.endpoint(q,'after')
        self.assertEqual(before.endpoint,q.segment.before)
        self.assertEqual(after.endpoint,q.segment.after)
        self.assertEqual(before.root,r); self.assertEqual(after.root,r)
        with self.assertRaises(ValueError): X.endpoint(q,'middle')

    def test_readiness_closes_ancestry_but_not_finite_master(self):
        d=X.readiness()
        self.assertTrue(d['correlated_COMPLETE_BRMM_left_inclusion_consumed'])
        self.assertTrue(d['primary_COMPLETE_BRMM_physical_condition_restriction_closed'])
        self.assertTrue(d['primary_COMPLETE_BRMM_requires_no_common_generator_representation'])
        self.assertTrue(d['bounded_primitive_maps_to_same_Live_origin_prefix_S'])
        self.assertTrue(d['bias_phi_driver_and_true_beta_hard_contracts_checked_per_segment'])
        self.assertTrue(d['source_cell_parent_child_and_primitive_continuity_checked'])
        self.assertFalse(d['qualified_async_endpoint_comes_from_admitted_transition'])
        self.assertTrue(d['qualified_async_endpoint_comes_from_checked_outer_transition'])
        self.assertTrue(d['sample_zero_checked_outer_endpoint_available_without_fake_transition'])
        self.assertFalse(d['sample_zero_full_source_membership_proved'])
        self.assertFalse(d['full_O601_membership_qualified_by_tokens_or_finite_checks'])
        self.assertTrue(d['raw_IMU_packet_bound_to_same_qualified_physical_predecessor'])
        self.assertTrue(d['persistent_gyro_and_accel_residual_history_tokens_required'])
        self.assertTrue(d['Racc_covariance_not_reinterpreted_as_hard_sensor_noise_bound'])
        self.assertFalse(d['quantitative_sensor_residual_ISS_envelope_attached'])
        self.assertFalse(d['finite_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(d['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(d['storage_search_allowed']); self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
