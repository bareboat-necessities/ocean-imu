"""Sample-zero magnetic/source interleave invariants for ALT.

These are finite necessary-source identities only.  They prove that a magnetic
call at the checked fresh-Live endpoint does not fabricate or consume a 5 ms
COMPLETE-BRMM transition, so transition ordinal 1 remains available to the first
source-bound IMU step.  They do not prove generator/QO membership or authorize
storage.
"""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_live_interleave as X
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import bias_families as BIAS
import test_finite_live_interleave as LIVEFIX
import test_finite_live_magnetic_word as MAG


def checked_origin(state):
    ref=state.live.live.mekf.reference
    root=X.SOURCE.certified_root(
        history_id=ref.history_id,
        generator_id='origin-regression-generator',
        live_origin=ref.live_origin,
        bias_family=ref.bias_family)
    return root,X.SOURCE.origin_endpoint(root,ref)


def mag_kwargs(state):
    kw=MAG.live_kwargs(state.magnetic)
    memory=state.magnetic.memory
    if memory.cfg.continuous_enabled:
        now=state.live.live.mekf.reference.time
        dt=now-memory.last_hi_time if memory.last_hi_time is not None and now>memory.last_hi_time else memory.cfg.sample_dt
        kw['hi_decay']=MAG.HI.Decay(dt,memory.cfg.continuous.memory,1)
    return kw


def family_phi(name):
    contract=next(c for c in BIAS.contracts() if c.name==name)
    return F.from_float(contract.phi_true.lo)


class Tests(unittest.TestCase):
    def test_sample_zero_mag_does_not_consume_first_physical_transition(self):
        s=LIVEFIX.root()
        root,origin=checked_origin(s)
        cont=X.SOURCE.begin(root)
        self.assertEqual(cont.next_ordinal,1)
        self.assertEqual(origin.endpoint,s.live.live.mekf.reference)

        mag=X.mag_step_source_qualified(s,origin,**mag_kwargs(s))
        self.assertEqual(cont.next_ordinal,1)
        self.assertEqual(mag.state.live.live.mekf.reference,origin.endpoint)
        self.assertEqual(mag.state.clock.calls,1)

        raw,segment,kw=LIVEFIX.imu_operands(mag.state)
        # LIVEFIX constructs a generic physical step.  The source-qualified
        # regression must use the actual fixed decay admitted by the selected
        # BIAS family rather than silently borrowing that generic fixture phi.
        segment=PHYS.PhysicalSegment(segment.before,segment.after,segment.J0,
            segment.J1,segment.J2,family_phi(root.bias_family),segment.bias_driver)
        witness=X.SOURCE.StepWitness(1,'root','cell-1','primitive-0','primitive-1')
        cont1=X.SOURCE.append(cont,witness=witness,segment=segment)
        sensor_root=X.SOURCE.SensorDisturbanceRoot(root,'gyro-residual-history','accel-residual-history')
        packet=X.SOURCE.qualify_raw_imu(cont1.steps[0],sensor_root,raw,'imu-packet-1')
        out=X.imu_step_source_qualified(mag.state,packet,**kw)

        self.assertEqual(cont1.steps[0].witness.ordinal,1)
        self.assertEqual(out.state.live.live.mekf.reference,segment.after)
        self.assertEqual(out.state.live.live.mekf.reference.time,origin.endpoint.time+F(1,200))

    def test_origin_object_is_necessary_outer_evidence_not_membership(self):
        s=LIVEFIX.root()
        _,origin=checked_origin(s)
        r=X.readiness()
        self.assertIsInstance(origin,X.SOURCE.QualifiedPhysicalOrigin)
        self.assertTrue(r['source_checked_async_magnetic_endpoint_entry_available'])
        self.assertFalse(r['sample_zero_full_source_membership_proved'])
        self.assertFalse(r['finite_magnetic_source_bound_to_same_COMPLETE_BRMM_history'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
