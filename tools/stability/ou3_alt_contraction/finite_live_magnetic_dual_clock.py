"""Dual-clock startup/Live magnetic composer for the ALT shipping word.

Shipping does not use one time coordinate in the outer magnetic path. Physical
source/inner-MEKF state advances on the IMU time base, while outer gating,
refinement and continuous hard-iron timestamps use ``SeaStateFusion_OU_III::t_``
(binary32). This module composes the existing magnetic algebra with the exact
canonical wrapper-clock descendant from ``finite_magnetic_wrapper_clock``.

The original finite_live_magnetic_word remains useful conditional finite-real
algebra. This stronger entry corrects its clock provenance from the FIRST
startup magnetic call, so accumulated sufficient statistics are never repaired
post hoc.  ``CertifiedStartupState`` is an inductive theorem object: it can only
be rooted at a pristine pre-magnetic startup state and only this module can
construct successors, each through the literal dual-clock startup_call below.
Binary32 exp/libm/Eigen correspondence is still open.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_live_magnetic_word as BASE
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as WCLOCK
from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as PREFIX

_CERT_TOKEN=object()


class CertifiedStartupState:
    """Inductive proof object for startup history executed only by this module.

    The constructor token deliberately prevents another theorem layer from
    wrapping an already-mutated BASE.StartupState and thereby laundering old
    single-clock ancestry into the dual-clock theorem path.
    """
    __slots__=('state','calls')
    def __init__(self,state,calls,*,_token=None):
        if _token is not _CERT_TOKEN:
            raise TypeError('dual-clock startup certificate is constructed only by certified entry/step')
        if not isinstance(state,BASE.StartupState):
            raise TypeError('BASE.StartupState required')
        if not isinstance(calls,int) or isinstance(calls,bool) or calls<0:
            raise ValueError('nonnegative certified startup call count required')
        self.state=state; self.calls=calls


def _pristine_startup(state:BASE.StartupState):
    """Require certificate creation before ANY magnetic startup mutation."""
    if not isinstance(state,BASE.StartupState):
        raise TypeError('persistent startup magnetic product state required')
    word=state.word; memory=state.memory
    if state.ready is not None:
        raise ValueError('dual-clock startup certificate must precede ready/reference capture')
    if (word.mag.tuner != BASE.TUNER.State() or word.mag.last_mag_time is not None
            or word.source_model_root is not None or word.source_history_id is not None):
        raise ValueError('dual-clock startup certificate cannot wrap an already-used magnetic word')
    if (memory.continuous != BASE.HI.State() or memory.last_hi_time is not None
            or memory.applied != BASE.HI.Applied()):
        raise ValueError('dual-clock startup certificate must precede continuous magnetic accumulation/application')
    return state


def certify_fresh_startup(state:BASE.StartupState):
    """Root the dual-clock induction at the literal pre-magnetic startup state."""
    return CertifiedStartupState(_pristine_startup(state),0,_token=_CERT_TOKEN)


@dataclass(frozen=True)
class CertifiedStartupResult:
    state: CertifiedStartupState
    event: BASE.StartupResult
    def __post_init__(self):
        if not isinstance(self.state,CertifiedStartupState) or not isinstance(self.event,BASE.StartupResult):
            raise TypeError('certified dual-clock startup successor and event required')
        if self.state.state != self.event.state:
            raise ValueError('certified startup successor detached from executed dual-clock event')


def _accumulate(memory, source, proxy, tilt, timestamp:WCLOCK.Timestamp, *,
                decay, eigen_success, spectrum):
    cfg=memory.cfg
    if not cfg.continuous_enabled or not proxy.initialized:
        if any(x is not None for x in (decay,eigen_success,spectrum)):
            raise ValueError('disabled/uninitialized continuous path consumes no witnesses')
        return memory,None
    if tilt is None:
        raise TypeError('continuous calibration needs the same private proxy tilt')
    dt=WCLOCK.shipping_elapsed(timestamp,memory.last_hi_time,fallback_dt=cfg.sample_dt)
    out=BASE.HI.update(memory.continuous,cfg.continuous,dt=dt,
                       rotation=BASE.CORE.rotation(tilt.q_tilt),raw_body=source.raw_body,
                       decay=decay,eigen_success=eigen_success,spectrum=spectrum)
    return replace(memory,continuous=out.state,last_hi_time=timestamp.wrapper_time),out


def startup_call(state:BASE.StartupState,physical,*,residual_body=None,packet_id=None,
                 begun=True,have_last_imu=True,proxy_q_norm=None,proxy_yaw_half=None,
                 hi_decay=None,hi_eigen_success=None,hi_spectrum=None,
                 mag_norm=None,mean_norm=None,horizontal_sqrt=None,ready_yaw_half=None):
    """Literal outer startup updateMag using binary32 wrapper time."""
    if not isinstance(state,BASE.StartupState):
        raise TypeError('persistent startup magnetic product state required')
    cfg=state.memory.cfg; ts=WCLOCK.at_physical_time(physical.time); wt=ts.wrapper_time
    if not isinstance(begun,bool) or not isinstance(have_last_imu,bool):
        raise TypeError('literal startup wrapper branches required')
    operands=(residual_body,packet_id,proxy_q_norm,proxy_yaw_half,hi_decay,
              hi_eigen_success,hi_spectrum,mag_norm,mean_norm,horizontal_sqrt,ready_yaw_half)
    if not begun or not cfg.gate.with_mag or wt < cfg.gate.mag_delay:
        if any(x is not None for x in operands):
            raise ValueError('outer-gated startup call consumes no source/arithmetic operands')
        return BASE.StartupResult(state,None,None,None)

    qualified=BASE._source(state.memory,physical,residual_body,packet_id)
    source=qualified.sample; proxy=state.word.mag.proxy
    need_tilt=cfg.continuous_enabled and proxy.initialized
    admission=BASE.GRAVITY.startup_mag_admission(
        state.word.gate,cfg.gravity,wrapper_time=wt,begun=begun,
        have_last_imu=have_last_imu,mag_ref_set=state.ready is not None)
    tilt=BASE._tilt(proxy,q_norm=proxy_q_norm,yaw_half=proxy_yaw_half) if need_tilt or admission.admitted else None
    if tilt is None and (proxy_q_norm is not None or proxy_yaw_half is not None):
        raise ValueError('unread private tilt consumes no quaternion witnesses')
    memory,continuous=_accumulate(state.memory,source,proxy,tilt,ts,
        decay=hi_decay,eigen_success=hi_eigen_success,spectrum=hi_spectrum)

    packet=PREFIX.Packet(wt,source.raw_body,source.packet_id)
    base=BASE.START.update_mag_call(
        state.word,cfg.gravity,cfg.tuner,packet,begun=begun,
        have_last_imu=have_last_imu,mag_ref_set=state.ready is not None,
        sample_dt=cfg.sample_dt,boat_q_norm=proxy_q_norm if admission.admitted else None,
        yaw_half=proxy_yaw_half if admission.admitted else None,
        mag_norm=mag_norm,mean_norm=mean_norm,horizontal_sqrt=horizontal_sqrt)
    root=source.model.model_root; history=source.physical.history_id
    if state.word.source_model_root is not None:
        if root!=state.word.source_model_root or history!=state.word.source_history_id:
            raise ValueError('startup magnetic source model/physical history restarted')
    word_state=BASE.START.State(base.state.gate,base.state.mag,
        state.word.source_model_root or root,state.word.source_history_id or history)
    out=BASE.START.MagResult(word_state,base.admission,base.magnetic,source,qualified)

    ready=state.ready
    due=(ready is None and out.magnetic is not None and out.magnetic.tuner_step.returned_ready
         and BASE.M.dot(out.state.mag.tuner.world_reference,out.state.mag.tuner.world_reference)>cfg.min_reference_norm**2)
    if due:
        ready=BASE.READY.transition(out.state.mag.tuner,
            BASE.READY.Config(cfg.sigma_internal,memory.producer_root),yaw_half=ready_yaw_half)
    elif ready_yaw_half is not None:
        raise ValueError('no provisional reference write consumes no gauge witness')
    return BASE.StartupResult(BASE.StartupState(out.state,memory,ready),out,continuous,qualified)


def certified_startup_call(cert:CertifiedStartupState, physical, **kwargs):
    """Advance the inductive startup certificate by one literal dual-clock call."""
    if not isinstance(cert,CertifiedStartupState):
        raise TypeError('certified dual-clock startup predecessor required')
    out=startup_call(cert.state,physical,**kwargs)
    nxt=CertifiedStartupState(out.state,cert.calls+1,_token=_CERT_TOKEN)
    return CertifiedStartupResult(nxt,out)


def live_call(state:BASE.LiveState,core,proxy,*,residual_body=None,packet_id=None,
              proxy_q_norm=None,proxy_yaw_half=None,
              hi_decay=None,hi_eigen_success=None,hi_spectrum=None,
              mag_norm=None,mean_norm=None,horizontal_sqrt=None,
              gauge_half=None,mekf_q_norm=None,mekf_yaw_half=None,
              apply_decay=None,apply_new_norm=None,apply_anchor_norm=None,
              ldlt=None,alpha=1,radius=None):
    """One Live magnetic call with outer and inner clocks kept distinct."""
    if not isinstance(state,BASE.LiveState) or not isinstance(core,BASE.CORE.State) or not isinstance(proxy,BASE.VERT.State):
        raise TypeError('persistent Live magnetic/core/private observer states required')
    cfg=state.memory.cfg; pt=core.reference.time; ts=WCLOCK.at_physical_time(pt); wt=ts.wrapper_time
    operands=(residual_body,packet_id,proxy_q_norm,proxy_yaw_half,hi_decay,
        hi_eigen_success,hi_spectrum,mag_norm,mean_norm,horizontal_sqrt,
        gauge_half,mekf_q_norm,mekf_yaw_half,apply_decay,apply_new_norm,apply_anchor_norm,ldlt,radius)
    if not cfg.gate.with_mag or wt < cfg.gate.mag_delay:
        if any(x is not None for x in operands) or alpha!=1:
            raise ValueError('outer-gated Live call consumes no source/arithmetic operands')
        return BASE.LiveResult(state,core,None,None,None,core,None,None,None,None)

    qualified=BASE._source(state.memory,core.reference,residual_body,packet_id)
    source=qualified.sample
    refining=cfg.refinement_enabled and not state.refinement_done and wt>=cfg.refinement_start
    need_tilt=refining or (cfg.continuous_enabled and proxy.initialized)
    tilt=BASE._tilt(proxy,q_norm=proxy_q_norm,yaw_half=proxy_yaw_half) if need_tilt else None
    if tilt is None and (proxy_q_norm is not None or proxy_yaw_half is not None):
        raise ValueError('unread private tilt consumes no quaternion witnesses')
    memory,continuous=_accumulate(state.memory,source,proxy,tilt,ts,
        decay=hi_decay,eigen_success=hi_eigen_success,spectrum=hi_spectrum)

    active,control=state.active,state.control
    tuner,last=state.tuner,state.last_mag_time
    started,done,done_time=state.refinement_started,state.refinement_done,state.refinement_time
    refinement=release=None
    if refining:
        if not started:
            tuner,last,started=BASE.TUNER.State(),None,True
        dt=WCLOCK.shipping_elapsed(ts,last,fallback_dt=cfg.sample_dt)
        last=wt
        tc=replace(cfg.tuner,min_window=cfg.refinement_window)
        corrected=tuple(source.raw_body[i]-memory.applied.total_bias[i] for i in range(3))
        refinement=BASE.TUNER.step(tuner,tc,q_tilt_bw=tilt.q_tilt,mag_body=corrected,dt=dt,
            q_norm=None if tuner.ready else BASE.TUNER.SqrtWitness(1,1),mag_norm=mag_norm,
            mean_norm=mean_norm,horizontal_sqrt=horizontal_sqrt)
        tuner=refinement.state
        if refinement.returned_ready and BASE.M.dot(tuner.world_reference,tuner.world_reference)>cfg.min_reference_norm**2:
            active=BASE.EVENTS.rewrite(active,BASE.EVENTS.WriteWitness('live_refinement',
                BASE.MAG.Model(tuner.world_reference,cfg.sigma_internal),memory.producer_root))
            core=BASE._overwrite_absolute_yaw(core,tuner.mean,gauge_half,mekf_q_norm,mekf_yaw_half)
            done,done_time=True,wt
            release=BASE.GATE.set_hold(control,core,cfg.gate,hold=False,live=True)
            control,core=release.state,release.filter_state
        elif any(x is not None for x in (gauge_half,mekf_q_norm,mekf_yaw_half)):
            raise ValueError('unfinished refinement consumes no MEKF yaw-write witnesses')
    elif any(x is not None for x in (mag_norm,mean_norm,horizontal_sqrt,gauge_half,mekf_q_norm,mekf_yaw_half)):
        raise ValueError('inactive refinement consumes no acquisition/yaw operands')

    after_refinement=core
    applied=BASE.HI.apply(memory.applied,memory.continuous,active.model.world_reference,
        time=wt,sample_dt=cfg.sample_dt,fraction=cfg.apply_fraction,slew_tau=cfg.slew_tau,
        enabled=cfg.continuous_enabled,refinement_enabled=cfg.refinement_enabled,
        refinement_done=done,proxy_min_norm=cfg.min_reference_norm,
        decay=apply_decay,new_norm=apply_new_norm,anchor_norm=apply_anchor_norm)
    memory=replace(memory,applied=applied.state)
    if applied.wrote_reference:
        active=BASE.EVENTS.rewrite(active,BASE.EVENTS.WriteWitness('continuous_hard_iron',
            BASE.MAG.Model(applied.reference,cfg.sigma_internal),memory.producer_root))
    corrected=tuple(source.raw_body[i]-memory.applied.total_bias[i] for i in range(3))
    field_difference=tuple(source.model.world_field[i]-active.model.world_reference[i] for i in range(3))
    rotated_difference=BASE.MAG.SENSOR.q_rotate(core.reference.q_world_to_body,field_difference)
    nu=tuple(rotated_difference[i]+source.model.hard_iron_body[i]
             -memory.applied.total_bias[i]+source.residual_body[i] for i in range(3))
    packet=BASE.MAG.Sample(core.reference,corrected,nu,active.model)
    measured=BASE.ASYNC.update_mag_call(BASE.ASYNC.State(core,control,active),cfg.gate,
        time=pt,live=True,sample=packet,ldlt=ldlt,alpha=alpha,radius=radius)
    nxt=BASE.LiveState(memory,tuner,last,active,measured.state.control,started,done,done_time)
    return BASE.LiveResult(nxt,measured.state.filter,qualified,continuous,refinement,
        after_refinement,release,applied,measured,nu)


def readiness():
    c=WCLOCK.readiness(); b=BASE.readiness()
    return {
      'existing_magnetic_event_order_reused':b['continuous_application_precedes_same_packet_measurement'],
      'startup_outer_delay_uses_exact_binary32_wrapper_clock':True,
      'startup_continuous_statistics_use_binary32_wrapper_elapsed_time':True,
      'startup_tuner_packet_clock_uses_binary32_wrapper_time':True,
      'fresh_startup_dual_clock_induction_object_available':True,
      'already_mutated_single_clock_startup_cannot_be_certified':True,
      'certified_startup_successor_only_from_dual_clock_call':True,
      'live_refinement_start_and_elapsed_use_binary32_wrapper_clock':True,
      'live_continuous_sample_and_apply_clocks_use_binary32_wrapper_clock':True,
      'inner_MEKF_magnetic_call_retains_physical_inner_time':True,
      'dual_clock_magnetic_word_composed':True,
      'canonical_prefix_wrapper_clock_arithmetic_closed':c['canonical_prefix_wrapper_clock_arithmetic_closed'],
      'binary32_exp_solver_roundoff_closed':False,
      'source_uniform_complete_magnetic_word_qualified':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
