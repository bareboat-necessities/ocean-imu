"""Compose the shipping H18 release guard with qualified Live mag timing.

This closes reachability of the *unlock guard*, not the entire H18->A21 word.
The shipping counter advances on every post-delay updateMag call.  Under
MAG-CALL-SCHEDULE-v1, the count reaches 250 within 10 s after Live. Call 250
need not satisfy the strict one-second guard. Continued, locally finite calls
force both predicates by max(249*gap, 1+gap) after the first call; on the
default gap this is 9.96 s, hence at most 10 s after Live.  An external acc-bias hold can still prevent immediate
A21 enable and therefore remains a separate hybrid/source obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as S

@dataclass(frozen=True)
class Result:
    mag_updates_applied:int
    first_mag_finite:bool
    elapsed_since_first_mag_upper:object
    live_to_unlock_upper:object
    accel_bias_lock_forced_clear:bool
    A21_enable_forced:bool


def reach(*,schedule:S.Schedule|None=None,external_hold:bool):
    r=S.release_reachability(schedule)
    S.require_release_reachable(r)
    # At some call by the corrected joint deadline: Live, count >= 250,
    # finite first-time and elapsed > 1 all hold. The actual call index may
    # exceed 250. mag_updates_applied records the guaranteed lower bound.
    clear=True
    return Result(r.unlock_count,True,r.elapsed_first_to_unlock_max,
                  r.live_to_unlock_max,clear,clear and not bool(external_hold))


def readiness():
    a=reach(external_hold=False);h=reach(external_hold=True)
    return {
      'qualified_mag_schedule_to_literal_unlock_guard_attached':True,
      'accel_bias_lock_forced_clear_within_10s_of_gauged_Live':a.accel_bias_lock_forced_clear,
      'unlock_reachability_independent_of_external_hold':h.accel_bias_lock_forced_clear,
      'A21_enable_closed_when_no_external_hold':a.A21_enable_forced,
      'A21_enable_closed_for_arbitrary_external_hold_history':False,
      'H18_A21_same_history_covariance_edge_composed_with_reachability':False,
      'ALT_LIVE_PASS':False,
    }
