"""Finite magnetic event identities; fixtures are not BRMM/storage evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_live_magnetic_word as X
from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as GAUGE
from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as PREFIX
from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as QUAL
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_continuous_mag_runtime as HI
import test_finite_startup_first_live_step as FIRST
import test_finite_continuous_mag_runtime as HF

PASS=MR.SafeLDLT(True,None,1,F(1,10**7))
REJECT=MR.SafeLDLT(False,False,1,F(1,10**7))


def zero_yaw(c=1):
    return X.TILT.YawHalfWitness(c,0,X.TILT.SqrtWitness(c*c,c),1,0)


def startup(*,refinement_start=90,refinement_window=0,continuous=True,hard_iron=(0,0,0)):
    bridge=FIRST.startup()
    cfg=X.Config(gate=X.GATE.Config(mag_delay=0),gravity=X.GRAVITY.Config(mag_delay=0),
                 tuner=X.TUNER.Config(min_samples=1,min_window=0),
                 refinement_start=refinement_start,refinement_window=refinement_window,
                 continuous_enabled=continuous)
    model=X.SOURCE.Model((30,0,0),hard_iron,'fixed-physical-model')
    proxy=bridge.frontend_before.tuner.vertical
    word=X.START.State(X.GRAVITY.State(gravity_good=2,aligned_branch=True),PREFIX.State(proxy))
    s=X.begin_startup(word,cfg,model,bridge.state.mekf.reference.history_id,'startup-reference')
    v=30+hard_iron[0]
    out=X.startup_call(s,bridge.state.mekf.reference,residual_body=(0,0,0),packet_id='startup',
        proxy_q_norm=X.TILT.SqrtWitness(1,1),proxy_yaw_half=zero_yaw(),
        hi_decay=HI.Decay(F(1,200),600,1) if continuous else None,
        mag_norm=X.TUNER.SqrtWitness(v*v,v),mean_norm=X.TUNER.SqrtWitness(v*v,v),
        horizontal_sqrt=GAUGE.HorizontalSqrt(v*v,v),ready_yaw_half=zero_yaw(v))
    return bridge,out.state


def live_kwargs(state,*,residual=(0,0,0),ldlt=REJECT):
    cfg=state.memory.cfg
    kw=dict(residual_body=residual,packet_id='live',ldlt=ldlt)
    if cfg.continuous_enabled:
        kw.update(proxy_q_norm=X.TILT.SqrtWitness(1,1),proxy_yaw_half=zero_yaw(),
                  hi_decay=HI.Decay(cfg.sample_dt,cfg.continuous.memory,1)
                  if cfg.continuous.memory>0 else None)
    return kw


def refine_kwargs(v=30):
    return dict(proxy_q_norm=X.TILT.SqrtWitness(1,1),proxy_yaw_half=zero_yaw(),
        mag_norm=X.TUNER.SqrtWitness(v*v,v),mean_norm=X.TUNER.SqrtWitness(v*v,v),
        horizontal_sqrt=GAUGE.HorizontalSqrt(v*v,v),gauge_half=zero_yaw(v),
        mekf_q_norm=X.TILT.SqrtWitness(1,1),mekf_yaw_half=zero_yaw())


class Tests(unittest.TestCase):
    def test_startup_packet_also_feeds_continuous_memory_before_Live(self):
        bridge,s=startup()
        self.assertEqual(s.memory.continuous.weight,1)
        self.assertEqual(s.memory.continuous.body_sum,(30,0,0))
        self.assertEqual(s.memory.last_hi_time,bridge.state.mekf.reference.time)
        self.assertEqual(s.ready.active_reference.model.world_reference,(30,0,0))
        live=X.enter_live(s)
        self.assertIs(live.memory,s.memory)
        self.assertEqual(live.control.updates,0); self.assertTrue(live.control.hold)

    def test_startup_gravity_rejection_does_not_disable_continuous_accumulation(self):
        bridge,s=startup()
        word=X.START.State(X.GRAVITY.State(),PREFIX.State(bridge.frontend_before.tuner.vertical))
        s=X.begin_startup(word,s.memory.cfg,s.memory.model,s.memory.history_id,s.memory.producer_root)
        out=X.startup_call(s,bridge.state.mekf.reference,residual_body=(0,0,0),packet_id='gated',
            proxy_q_norm=X.TILT.SqrtWitness(1,1),proxy_yaw_half=zero_yaw(),
            hi_decay=HI.Decay(F(1,200),600,1))
        self.assertFalse(out.startup.admission.admitted)
        self.assertEqual(out.state.memory.continuous.weight,1)
        self.assertEqual(out.state.word.mag.tuner.accumulator.accepted_count,0)

    def test_continuous_history_cannot_be_reset_at_startup_handoff(self):
        _,s=startup()
        with self.assertRaisesRegex(ValueError,'restarted'):
            X.begin_startup(s.word,s.memory.cfg,s.memory.model,s.memory.history_id,s.memory.producer_root)

    def test_same_sample_subtraction_and_reference_discrepancy_are_not_free_noise(self):
        bridge,s=startup(hard_iron=(2,0,0)); live=X.enter_live(s)
        out=X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,
                        **live_kwargs(live,residual=(1,1,0)))
        self.assertEqual(out.measurement.magnetic.sample.raw_body,(33,1,0))
        self.assertEqual(out.effective_residual,(1,1,0))  # B_true=30, learned=32, physical HI=2
        self.assertEqual(out.qualification.sample.raw_body,(33,1,0))
        self.assertEqual(out.state.memory.continuous.body_sum,(65,1,0))
        self.assertFalse(out.measurement.measurement_accepted)
        self.assertEqual(out.state.control.updates,1)

    def test_fixed_source_model_cannot_be_supplied_per_call(self):
        bridge,s=startup(); live=X.enter_live(s)
        with self.assertRaises(TypeError):
            X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,model=s.memory.model)
        bad=replace(bridge.state.mekf,reference=replace(bridge.state.mekf.reference,history_id='other'))
        with self.assertRaisesRegex(ValueError,'physical history'):
            X.live_call(live,bad,bridge.state.tuner.vertical,**live_kwargs(live))

    def test_out_of_source_envelope_fails_before_magnetic_word(self):
        bridge,s=startup(); live=X.enter_live(s)
        with self.assertRaisesRegex(ValueError,'mag residual exceeds'):
            X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,
                        **live_kwargs(live,residual=(3,0,0)))
        self.assertEqual(live.memory.continuous.weight,1)

    def test_squared_source_qualification_accepts_irrational_norm_without_guess(self):
        p=X.SOURCE.PhysicalEndpoint(0,(1,0,0,0),'h')
        model=X.SOURCE.Model((30,1,0),(1,1,0),'m')
        packet=X.SOURCE.make_sample(p,model,(1,1,0),'p')
        q=QUAL.qualify_squared(packet)
        self.assertTrue(q.qualified)
        bad=X.SOURCE.make_sample(p,model,(2,1,0),'bad')
        with self.assertRaisesRegex(ValueError,'exceeds'):
            QUAL.qualify_squared(bad)

    def test_refinement_resets_acquisition_and_writes_before_same_sample_measurement(self):
        bridge,s=startup(refinement_start=0); live=X.enter_live(s)
        kw=live_kwargs(live,ldlt=PASS); kw.update(refine_kwargs())
        out=X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,**kw)
        self.assertEqual(out.refinement.state.accumulator.accepted_count,1)
        self.assertTrue(out.state.refinement_done); self.assertFalse(out.state.control.hold)
        self.assertEqual(out.state.active.generation,1)
        self.assertEqual(out.measurement.state.magnetic_active,out.state.active)
        self.assertEqual(out.refinement_filter.covariance,bridge.state.mekf.covariance)
        self.assertTrue(out.measurement.measurement_accepted)
        self.assertEqual(out.state.control.updates,1)
        self.assertEqual(out.state.memory.continuous.weight,2)
        self.assertFalse(out.continuous_apply.wrote_reference)  # no fitted estimate yet

    def test_refinement_reject_keeps_clock_and_prior_active_reference(self):
        bridge,s=startup(refinement_start=0,refinement_window=F(1,100)); live=X.enter_live(s)
        kw=live_kwargs(live);kw.update(mag_norm=X.TUNER.SqrtWitness(900,30))
        first=X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,**kw)
        self.assertFalse(first.state.refinement_done)
        self.assertEqual(first.state.last_mag_time,0)
        self.assertEqual(first.state.active,live.active)
        # Same timestamp uses fallback .005 again and reaches .010 accepted window.
        kw=live_kwargs(first.state);kw.update(refine_kwargs())
        second=X.live_call(first.state,first.filter,bridge.state.tuner.vertical,**kw)
        self.assertTrue(second.state.refinement_done)
        self.assertEqual(second.refinement.state.accumulator.accepted_window,F(1,100))

    def test_detached_refinement_yaw_witness_fails(self):
        bridge,s=startup(refinement_start=0);live=X.enter_live(s)
        kw=live_kwargs(live);kw.update(refine_kwargs());kw['gauge_half']=zero_yaw(31)
        with self.assertRaisesRegex(ValueError,'same accepted mean'):
            X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,**kw)

    def test_nondegenerate_yaw_write_keeps_all_nonattitude_coordinates_and_full_covariance(self):
        bridge,_=startup()
        core=bridge.state.mekf
        # Mean direction (7,24) has half-angle (4/5,3/5). No trace fitting.
        gauge=X.TILT.YawHalfWitness(7,24,X.TILT.SqrtWitness(625,25),F(4,5),F(3,5))
        new=X._overwrite_absolute_yaw(core,(7,24,0),gauge,X.TILT.SqrtWitness(1,1),zero_yaw())
        self.assertEqual(new.q_hat,(F(4,5),0,0,F(3,5)))
        self.assertEqual(new.z[:3],(0,0,F(-3,2)))
        self.assertEqual(new.z[3:],core.z[3:]);self.assertEqual(new.covariance,core.covariance)
        self.assertIs(new.reference,core.reference)

    def test_refinement_releases_previously_unlocked_hold_before_measurement(self):
        bridge,s=startup(refinement_start=0); live=X.enter_live(s)
        # Conditional control-edge fixture: not a source-uniform startup claim.
        live=replace(live,control=X.GATE.State(250,0,False,True))
        kw=live_kwargs(live);kw.update(refine_kwargs())
        out=X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,**kw)
        self.assertEqual(out.refinement_filter.mode,'A')
        self.assertTrue(out.hold_release.enabled_bias_now)
        self.assertEqual(out.filter.mode,'A')
        self.assertFalse(out.measurement.measurement_accepted)

    def test_refinement_and_continuous_rewrite_can_both_happen_on_one_call(self):
        bridge,s=startup(refinement_start=0,hard_iron=(2,0,0));live=X.enter_live(s)
        stats,hcfg=HF.fit(True)
        # Local event identity at conditional carried statistics; not promoted
        # to a startup-reachable or BRMM-qualified initial state.
        cfg=replace(live.memory.cfg,continuous=hcfg,slew_tau=0)
        live=replace(live,memory=replace(live.memory,cfg=cfg,continuous=stats))
        kw=live_kwargs(live);kw.update(refine_kwargs(32))
        kw.update(apply_new_norm=HI.HorizontalNorm((F(91,3),0),F(91,3)),
                  apply_anchor_norm=HI.HorizontalNorm((F(92,3),0),F(92,3)))
        out=X.live_call(live,bridge.state.mekf,bridge.state.tuner.vertical,**kw)
        self.assertEqual(out.state.active.generation,2)
        self.assertEqual(out.state.memory.applied.applied,(1,0,0))
        self.assertEqual(out.state.active.model.world_reference,(F(95,3),0,0))
        self.assertEqual(out.measurement.magnetic.sample.raw_body,(31,0,0))
        self.assertEqual(out.effective_residual,(F(-2,3),0,0))
        self.assertEqual(out.measurement.state.magnetic_active,out.state.active)
        self.assertEqual(out.state.active.model.sigma_internal,live.active.model.sigma_internal)

    def test_readiness_does_not_promote_conditional_component_words(self):
        r=X.readiness()
        for key in ('complete_word_finite_identity','ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS',
                    'deployment_finite_precision_closed','source_uniform_complete_magnetic_word_qualified'):
            self.assertFalse(r[key])


if __name__=='__main__':unittest.main()
