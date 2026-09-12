"""Deterministic asynchronous magnetometer call schedule for ALT.

Shipping increments ``mag_updates_applied_`` on every post-delay updateMag call
immediately after invoking the MEKF magnetometer update; the count is not gated
on innovation acceptance.  ALT therefore admits a deployment call schedule,
separate from magnetic-value/source qualification. Calls are locally finite
and continue over unbounded physical time (not a Zeno list):

    after magnetically gauged Live, first updateMag call <= 0.04 s,
    and every subsequent updateMag call gap <= 0.04 s.

This is a deployment/source timing assumption, not a BMM150 electrical guarantee.
It is only 25 Hz, intentionally weaker than the target sensor capability and the
same cadence used by the shipping finite-word host regression.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

ASSUMPTION_ID='MAG-CALL-SCHEDULE-v1'
DEFAULT_MAX_GAP=F(1,25)
DEFAULT_UNLOCK_COUNT=250
DEFAULT_STRICT_GUARD=F(1)

@dataclass(frozen=True)
class Schedule:
    first_after_live_max:F=DEFAULT_MAX_GAP
    gap_max:F=DEFAULT_MAX_GAP
    assumption_id:str=ASSUMPTION_ID
    def __post_init__(self):
        a,b=F(self.first_after_live_max),F(self.gap_max)
        if a<=0 or b<=0: raise ValueError('positive deterministic magnetometer timing bounds required')
        if self.assumption_id!=ASSUMPTION_ID: raise ValueError('wrong magnetometer call-schedule assumption')
        object.__setattr__(self,'first_after_live_max',a);object.__setattr__(self,'gap_max',b)

@dataclass(frozen=True)
class Reachability:
    unlock_count:int
    first_call_after_live_max:F
    elapsed_first_to_unlock_max:F
    live_to_unlock_max:F
    strict_one_second_guard_satisfied:bool


def default_schedule(): return Schedule()


def release_reachability(schedule:Schedule|None=None,unlock_count=DEFAULT_UNLOCK_COUNT):
    s=default_schedule() if schedule is None else schedule
    if not isinstance(s,Schedule): raise TypeError('qualified magnetometer call schedule required')
    if not isinstance(unlock_count,int) or isinstance(unlock_count,bool) or unlock_count<1:
        raise ValueError('positive integer unlock count required')
    n=unlock_count
    # An upper inter-call bound supplies NO lower bound on t_(n-1)-t_0.
    # Count n may arrive before the strict guard. In that case the first call
    # strictly after t_0+1 already has count >= n and occurs within one gap.
    # Existence uses continued, time-unbounded, locally finite call coverage.
    count_deadline=(n-1)*s.gap_max
    guard_deadline=DEFAULT_STRICT_GUARD+s.gap_max
    elapsed=max(count_deadline,guard_deadline)
    total=s.first_after_live_max+elapsed
    strict=True
    return Reachability(n,s.first_after_live_max,elapsed,total,strict)


def require_release_reachable(r:Reachability):
    if not isinstance(r,Reachability) or not r.strict_one_second_guard_satisfied:
        raise ValueError('call schedule does not force the strict one-second unlock guard')
    return True


@dataclass(frozen=True)
class Clock:
    live_time:F
    last_time:F|None=None
    calls:int=0
    def __post_init__(self):
        t=F(self.live_time)
        if t<0 or not isinstance(self.calls,int) or isinstance(self.calls,bool) or self.calls<0:
            raise ValueError('valid persistent Live magnetic clock required')
        object.__setattr__(self,'live_time',t)
        if self.last_time is not None:
            last=F(self.last_time)
            if last<t: raise ValueError('magnetic call cannot precede Live')
            object.__setattr__(self,'last_time',last)
        if (self.calls==0)!=(self.last_time is None):
            raise ValueError('magnetic count and last-call clock must agree')


def check_prefix(clock:Clock,schedule:Schedule,*,time):
    """Validate finite-prefix timing only; never infer infinite coverage."""
    if not isinstance(clock,Clock) or not isinstance(schedule,Schedule):
        raise TypeError('persistent magnetic clock and schedule required')
    t=F(time)
    start=clock.live_time if clock.last_time is None else clock.last_time
    gap=schedule.first_after_live_max if clock.last_time is None else schedule.gap_max
    if t<start: raise ValueError('magnetic physical clock moved backwards')
    if t>start+gap: raise ValueError('declared magnetic call-schedule deadline missed')
    return True


def record_call(clock:Clock,schedule:Schedule,*,time):
    check_prefix(clock,schedule,time=time)
    return Clock(clock.live_time,F(time),clock.calls+1)


def readiness():
    r=release_reachability()
    return {
      'named_deterministic_mag_call_schedule_declared':True,
      'first_Live_mag_call_max_gap_s':r.first_call_after_live_max,
      'subsequent_mag_call_max_gap_s':DEFAULT_MAX_GAP,
      'default_250_count_elapsed_from_first_s':r.elapsed_first_to_unlock_max,
      'default_250_count_reached_within_Live_s':r.live_to_unlock_max,
      'strict_one_second_guard_forced':r.strict_one_second_guard_satisfied,
      'count_threshold_call_itself_need_not_satisfy_time_guard':True,
      'locally_finite_time_unbounded_call_coverage_required':True,
      'innovation_acceptance_not_used_for_shipping_counter':True,
      'external_acc_bias_hold_excluded_here':False,
      'H18_A21_complete_word_edge_attached':False,
      'ALT_LIVE_PASS':False,
    }
