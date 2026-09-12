"""Source-owning Live product regressions; not stability evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as WCLOCK
from tools.stability.ou3_alt_contraction import proof_plan as PLAN
import test_finite_live_interleave as BASE

BIAS0=next(c for c in BIAS.contracts() if c.name=='BIAS0')
PHI=F.from_float(BIAS0.phi_true.lo)


def root_state():
    live=BASE.root(); ref=live.live.live.mekf.reference
    root=SOURCE.SourceRoot(ref.history_id,'generator',ref.live_origin,'BIAS0',BIAS0.parameter_token)
    sensors=SOURCE.SensorDisturbanceRoot(root,'gyro-history','accel-history')
    _,_,kw=BASE.imu_operands(live)
    runtime=X.RuntimeConfig.from_step_kwargs(kw)
    return X.from_live(live,root,sensors,runtime)


def next_imu_operands(state):
    raw,seg,kw=BASE.imu_operands(state.live)
    qualified_seg=PHYS.PhysicalSegment(seg.before,seg.after,seg.J0,
        seg.J1,seg.J2,PHI,seg.bias_driver)
    k=state.source.next_ordinal
    parent='root' if k==1 else state.source.steps[-1].witness.source_cell_id
    pin='p0' if k==1 else state.source.steps[-1].witness.primitive_out_id
    witness=SOURCE.StepWitness(k,parent,f'c{k}',pin,f'p{k}')
    dynamic=X.dynamic_kwargs(kw); dynamic.pop('dt',None)
    return witness,qualified_seg,raw,dynamic


def imu(state):
    witness,segment,raw,kw=next_imu_operands(state)
    return X.imu_step(state,witness=witness,segment=segment,raw=raw,
                      packet_id=f'imu-{witness.ordinal}',**kw)


def mag_kwargs(state):
    kw=BASE.MAG.live_kwargs(state.live.magnetic)
    memory=state.live.magnetic.memory; physical=state.live.live.live.mekf.reference.time
    if memory.cfg.continuous_enabled:
        ts=WCLOCK.at_physical_time(physical)
        dt=WCLOCK.shipping_elapsed(ts,memory.last_hi_time,fallback_dt=memory.cfg.sample_dt)
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
        self.assertEqual(second.state.bias_history_id,s.bias_history_id)
        self.assertIs(second.state.runtime,s.runtime)

    def test_post_first_IMU_mag_event_cannot_advance_or_restart_source(self):
        s=imu(root_state()).state
        before=s.source
        event=X.mag_step(s,**mag_kwargs(s))
        self.assertIs(event.state.source,before)
        self.assertIs(event.state.runtime,s.runtime)
        self.assertEqual(event.state.live.live.live.mekf.reference,before.steps[-1].segment.after)
        self.assertEqual(event.state.source.root,before.root)
        self.assertEqual(event.state.bias_history_id,s.bias_history_id)

    def test_sample_zero_mag_uses_actual_fresh_origin_without_advancing_source(self):
        s=root_state(); before=s.source
        event=X.mag_step(s,**mag_kwargs(s))
        self.assertIs(event.state.source,before)
        self.assertEqual(len(event.state.source.steps),0)
        self.assertEqual(event.state.source.next_ordinal,1)
        self.assertEqual(event.state.live.live.live.mekf.reference.time,s.source.root.live_origin)
        first=imu(event.state)
        self.assertEqual(first.state.source.steps[0].witness.ordinal,1)

    def test_sample_zero_mag_rejects_fresh_reference_outside_physical_outer_cap(self):
        s=root_state(); core=s.live.live.live.mekf
        bad=replace(core,reference=replace(core.reference,acceleration=(100,0,0)))
        live=replace(s.live,live=replace(s.live.live,live=replace(s.live.live.live,mekf=bad)))
        q=replace(s,live=live)
        with self.assertRaisesRegex(ValueError,'acceleration vector cap'):
            X.mag_step(q,**mag_kwargs(q))
        self.assertEqual(q.live.magnetic.control.updates,0)

    def test_physical_source_constraint_runs_before_shipping_event(self):
        s=root_state(); witness,q,raw,kw=next_imu_operands(s)
        q=replace(q,after=replace(q.after,acceleration=(100,0,0)))
        with patch.object(X.LIVE,'imu_step_source_qualified') as event:
            with self.assertRaisesRegex(ValueError,'acceleration vector cap'):
                X.imu_step(s,witness=witness,segment=q,raw=raw,packet_id='invalid-source',**kw)
            event.assert_not_called()
        self.assertEqual(len(s.source.steps),0)

    def test_wrong_ordinal_fails_before_shipping_event(self):
        s=root_state(); _,segment,raw,kw=next_imu_operands(s)
        bad=SOURCE.StepWitness(2,'root','c2','p0','p2')
        with self.assertRaisesRegex(ValueError,'next source ordinal'):
            X.imu_step(s,witness=bad,segment=segment,raw=raw,packet_id='bad',**kw)
        self.assertEqual(len(s.source.steps),0)
        self.assertEqual(s.live.live.live.mekf.reference.time,s.source.root.live_origin)

    def test_static_runtime_config_and_dt_cannot_be_overridden_per_event(self):
        s=root_state(); witness,segment,raw,kw=next_imu_operands(s)
        with self.assertRaisesRegex(TypeError,'cannot override carried runtime config'):
            X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='bad',
                       commit_cfg=s.runtime.commit_cfg,**kw)
        with self.assertRaisesRegex(TypeError,'dt is owned'):
            X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='bad',dt=segment.h,**kw)

    def test_hold_event_preserves_source_continuation(self):
        s=imu(root_state()).state
        out=X.set_hold(s,hold=False)
        self.assertIs(out.state.source,s.source)
        self.assertIs(out.state.sensor_root,s.sensor_root)
        self.assertIs(out.state.runtime,s.runtime)
        self.assertEqual(out.state.bias_history_id,s.bias_history_id)
        self.assertEqual(out.state.live.live.live.mekf.reference,s.source.steps[-1].segment.after)

    def test_detached_sensor_root_is_rejected_by_product_state(self):
        s=root_state(); other=SOURCE.SourceRoot(s.source.root.history_id,'other-generator',
            s.source.root.live_origin,'BIAS0',BIAS0.parameter_token)
        sensors=SOURCE.SensorDisturbanceRoot(other,'g','a')
        with self.assertRaisesRegex(ValueError,'sensor histories detached'):
            X.State(s.live,s.source,sensors,s.bias_history_id,s.runtime)

    def test_same_family_but_restarted_bias_history_is_rejected(self):
        s=root_state(); core=s.live.live.live.mekf
        changed=replace(core,reference=replace(core.reference,bias_root='new-bias-history'))
        live=replace(s.live,live=replace(s.live.live,live=replace(s.live.live.live,mekf=changed)))
        with self.assertRaisesRegex(ValueError,'bias history restarted'):
            X.State(live,s.source,s.sensor_root,s.bias_history_id,s.runtime)

    def test_finite_master_status_closes_only_physical_forcing(self):
        status=X.finite_storage_status()
        self.assertEqual(status['map_representation'],'finite_physical_descriptor')
        self.assertTrue(status['physical_reference_forcing_retained'])
        self.assertTrue(status['zero_wind_heel_scope_enforced'])
        self.assertFalse(status['finite_error_identity_for_every_event'])
        self.assertFalse(status['all_coefficient_product_graphs_retained'])
        self.assertFalse(status['all_configured_branches_bound_to_finite_graph'])
        with self.assertRaisesRegex(RuntimeError,'finite-state storage blocked'):
            PLAN.assert_finite_storage_master(status)

    def test_readiness_advances_structure_without_promoting_master(self):
        r=X.readiness()
        self.assertTrue(r['source_continuation_is_part_of_theorem_product_state'])
        self.assertTrue(r['every_IMU_event_appends_exactly_next_source_ordinal'])
        self.assertTrue(r['analytic_BIAS_family_token_and_concrete_bias_history_both_persist'])
        self.assertTrue(r['persistent_static_runtime_configuration_carried_in_product_state'])
        self.assertTrue(r['theorem_IMU_event_cannot_override_static_runtime_configuration'])
        self.assertTrue(r['theorem_IMU_dt_owned_by_qualified_physical_segment'])
        self.assertTrue(r['magnetic_and_hold_events_preserve_current_source_endpoint'])
        self.assertTrue(r['physical_reference_forcing_retained_in_finite_master_status'])
        self.assertTrue(r['sample_zero_startup_to_checked_outer_endpoint_bridge_closed'])
        self.assertTrue(r['sample_zero_magnetic_entry_uses_checked_fresh_physical_origin'])
        self.assertTrue(r['sample_zero_magnetic_entry_advances_no_source_ordinal'])
        self.assertFalse(r['sample_zero_full_source_membership_proved'])
        self.assertFalse(r['sample_zero_startup_to_COMPLETE_BRMM_endpoint_bridge_closed'])
        self.assertFalse(r['finite_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
