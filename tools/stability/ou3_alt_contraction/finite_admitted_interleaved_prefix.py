"""Finite interleaving induction over the strongest admitted-source ALT edges.

The product persistently carries the canonical t_L restriction of the same
admitted COMPLETE-BRMM history.  That origin is checked against the fresh MEKF
physical Reference at ``begin`` and transition 1 must start from it.  Thereafter
all IMU/MAG/HOLD successors preserve the origin and both admitted histories.

IMU uses the canonical admitted BRMM+BIAS strong edge with source-bound runtime
coefficients and retained forcing; MAG uses the admitted magnetic edge and
retained correlated forcing; HOLD changes no physical-source ordinal.

All deployment arithmetic witnesses and startup reachability remain open, so a
600-transition trace is structurally representable but not a certified
source-uniform shipping word. No storage search is authorized.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as LIVE
from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_source_mag_word as MAG
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ABRMM
from tools.stability.ou3_alt_contraction import finite_admitted_startup_origin as START
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE


@dataclass(frozen=True)
class EventRecord:
    kind: str
    source_steps_before: int
    source_steps_after: int
    def __post_init__(self):
        if self.kind not in ('imu','mag','hold'):
            raise ValueError('unknown interleaved event kind')
        expected=self.source_steps_before+(1 if self.kind=='imu' else 0)
        if self.source_steps_after != expected:
            raise ValueError('event/source ordinal accounting mismatch')


@dataclass(frozen=True)
class State:
    live: LIVE.State
    origin: ABRMM.RestrictedOrigin
    events: tuple[EventRecord,...]=()
    def __post_init__(self):
        if not isinstance(self.live,LIVE.State) or not isinstance(self.origin,ABRMM.RestrictedOrigin):
            raise TypeError('joint admitted Live state and canonical t_L origin required')
        object.__setattr__(self,'events',tuple(self.events))
        if self.origin.history != self.live.admitted_history:
            raise ValueError('interleaved origin detached from carried admitted history')
        root=self.live.live_word.source.root
        ABRMM.qualify_origin(root,self.origin)
        n=len(self.live.live_word.source.steps)
        if n>SOURCE.TRANSITIONS: raise ValueError('more than canonical 600 source transitions')
        if n:
            if self.live.live_word.source.steps[0].segment.before != self.origin.endpoint:
                raise ValueError('interleaved source chain did not start at admitted t_L origin')
        else:
            ref=self.live.live_word.live.live.live.mekf.reference
            if ref != self.origin.endpoint:
                raise ValueError('fresh interleaved Reference is not admitted t_L origin')
        if self.events:
            if self.events[-1].source_steps_after != n:
                raise ValueError('event ledger detached from carried source continuation')
            if sum(e.kind=='imu' for e in self.events) != n:
                raise ValueError('IMU event count detached from source transition count')

    @property
    def imu_steps(self): return len(self.live.live_word.source.steps)
    @property
    def complete_source_horizon(self): return self.imu_steps==SOURCE.TRANSITIONS


def begin(live:LIVE.State, origin:ABRMM.RestrictedOrigin):
    # Reuse the dedicated startup binder so the fresh joint24 physical Reference
    # and canonical admitted-history origin have one identity check.
    START.bind(live,origin)
    return State(live,origin,())


def imu_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    restricted=kwargs.get('restricted')
    if not isinstance(restricted,ABRMM.RestrictedSegment):
        raise TypeError('interleaved IMU requires admitted BRMM restriction')
    if state.imu_steps==0 and restricted.segment.before != state.origin.endpoint:
        raise ValueError('first interleaved IMU restriction does not start at admitted t_L origin')
    before=state.imu_steps
    out=IMU.imu_step(state.live,**kwargs)
    record=EventRecord('imu',before,len(out.state.live_word.source.steps))
    return State(out.state,state.origin,state.events+(record,)),out


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    if 'origin' in kwargs:
        raise TypeError('interleaved prefix owns the admitted origin')
    before=state.imu_steps
    mkw=dict(kwargs)
    if before==0: mkw['origin']=state.origin
    out=MAG.mag_step(state.live,**mkw)
    record=EventRecord('mag',before,len(out.state.live_word.source.steps))
    return State(out.state,state.origin,state.events+(record,)),out


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    before=state.imu_steps
    out=LIVE.set_hold(state.live,hold=hold)
    record=EventRecord('hold',before,len(out.state.live_word.source.steps))
    return State(out.state,state.origin,state.events+(record,)),out


def readiness():
    i=IMU.readiness(); m=MAG.readiness(); s=START.readiness()
    return {
      'finite_successor_induction_over_IMU_MAG_HOLD_closed':True,
      'canonical_admitted_tL_origin_persistent_in_interleaved_product':True,
      'startup_sample_zero_identity_consumed_before_interleaving':s['startup_sample_zero_equal_to_admitted_history_restriction_proved'],
      'first_IMU_forced_to_start_at_same_admitted_tL_origin':True,
      'only_IMU_advances_admitted_source_ordinal':True,
      'IMU_successor_uses_strong_admitted_source_edge':i['physical_moments_sensor_forcing_and_model_roots_share_one_event'],
      'MAG_successor_retains_same_event_correlated_forcing':m['same_event_correlated_magnetic_forcing_retained'],
      'sample_zero_MAG_uses_persistent_admitted_origin':True,
      'both_admitted_histories_and_origin_persist_through_all_event_types':True,
      'canonical_600_transition_source_horizon_representable':True,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'sensor_disturbance_admissibility_attached':False,
      'startup_reachability_to_admitted_fresh_state_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
