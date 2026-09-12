"""Source-owning Live product regressions; not stability evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
import test_finite_live_interleave as BASE

BIAS0=next(c for c in BIAS.contracts() if c.name=='BIAS0')
PHI=F(str(BIAS0.phi_true.lo))


def root_state():
    live=BASE.root(); ref=live.live.live.mekf.reference
    root=SOURCE.SourceRoot(ref.history_id,'generator',ref.live_origin,'BIAS0',BIAS0.parameter_token)
    sensors=SOURCE.SensorDisturbanceRoot(root,'gyro-history','accel-history')
    return X.from_live(live,root,sensors)


def next_imu_operands(state):
    raw,seg,kw=BASE.imu_operands(state.live)
    # The component fixture uses phi=1 because beta=0 makes that algebraically
    # valid. Source admission must instead use the selected BIAS-family phi.
    qualified_seg=PHYS.PhysicalSegment(seg.before,seg.after,seg.delta_theta,
        seg.delta_velocity,seg.delta_position,PHI,seg.bias_driver)
    k=state.source.next_ordinal
    parent='root' if k==1 else state.source.steps[-1].witness.source_cell_id
    pin='p0' if k==1 else state.source.steps[-1].witness.primitive_out_id
    witness=SOURCE.StepWitness(k,parent,f'c{k}',pin,f'p{k}')
    return witness,qualified_seg,raw,kw


def imu(state):
    witness,segment,raw,kw=next_imu_operands(state)
    return X.imu_step(state,witness=witness,segment=segment,raw=raw,
                      packet_id=f'imu-{witness.ordinal}',**kw)


def mag_kwargs(state):
    kw=BASE.MAG.live_kwargs(state.live.magnetic)
    memory=state.live.magnetic.memory; now=state.live.live.live.mekf.reference.time
    if memory.cfg.continuous_enabled:
        dt=now-memory.last_hi_time if memory.last_hi_time is not None and now>memory.last_hi_time else memory.cfg.sample_dt
        kw['hi_decay']=BASE.MAG.HI.Decay(dt,memory.cfg.continuous.memory,1)
    return kw


class Tests(unittest.TestCase):
    def test_source_chain_advances_with_successive_IMU_events(self):
        s=root_state(); first=imu(s)
        self.assertEqual(len(first.state.source.steps),1)
        self.assertEqual(first.state.source.steps[-1].segment.after,
                         first.state.live.live.live.mekf.reference)
        second=imu(first.state)
        self.assertEqual(len(second.state.source.steps),2)
        self.assertEqual(second.state.source.steps[1].segment.before,
                         second.state.source.steps[0].segment.after)
        self.assertEqual(second.state.source.steps[1].witness.parent_source_cell_id,'c1')
        self.assertEqual(second.state.source.steps[1].witness.primitive_in_id,'p1')
        self.assertEqual(second.state.live.live.live.mekf.reference.time,F(1,100))

    def test_post_first_IMU_mag_event_cannot_advance_or_restart_source(self):
        s=imu(root_state()).state
        before=s.source
        event=X.mag_step(s,**mag_kwargs(s))
        self.assertIs(event.state.source,before)
        self.assertEqual(event.state.live.live.live.mekf.reference,before.steps[-1].segment.after)
        self.assertEqual(event.state.source.root,before.root)

    def test_sample_zero_mag_bridge_stays_explicitly_open(self):
        s=root_state()
        with self.assertRaisesRegex(NotImplementedError,'sample-zero'):
            X.mag_step(s,**mag_kwargs(s))

    def test_wrong_ordinal_fails_before_shipping_event(self):
        s=root_state(); _,segment,raw,kw=next_imu_operands(s)
        bad=SOURCE.StepWitness(2,'root','c2','p0','p2')
        with self.assertRaisesRegex(ValueError,'next source ordinal'):
            X.imu_step(s,witness=bad,segment=segment,raw=raw,packet_id='bad',**kw)
        self.assertEqual(len(s.source.steps),0)
        self.assertEqual(s.live.live.live.mekf.reference.time,s.source.root.live_origin)

    def test_hold_event_preserves_source_continuation(self):
        s=imu(root_state()).state
        out=X.set_hold(s,hold=False)
        self.assertIs(out.state.source,s.source)
        self.assertIs(out.state.sensor_root,s.sensor_root)
        self.assertEqual(out.state.live.live.live.mekf.reference,s.source.steps[-1].segment.after)

    def test_detached_sensor_root_is_rejected_by_product_state(self):
        s=root_state(); other=SOURCE.SourceRoot(s.source.root.history_id,'other-generator',
            s.source.root.live_origin,'BIAS0',BIAS0.parameter_token)
        sensors=SOURCE.SensorDisturbanceRoot(other,'g','a')
        with self.assertRaisesRegex(ValueError,'sensor histories detached'):
            X.State(s.live,s.source,sensors)

    def test_readiness_advances_structure_without_promoting_master(self):
        r=X.readiness()
        self.assertTrue(r['source_continuation_is_part_of_theorem_product_state'])
        self.assertTrue(r['every_IMU_event_appends_exactly_next_source_ordinal'])
        self.assertTrue(r['magnetic_and_hold_events_preserve_current_source_endpoint'])
        self.assertFalse(r['sample_zero_startup_to_COMPLETE_BRMM_endpoint_bridge_closed'])
        self.assertFalse(r['finite_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
