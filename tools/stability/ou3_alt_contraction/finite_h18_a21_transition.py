"""Compose qualified mag-call reachability with the literal H18/A21 hybrid edge.

The shipping theorem must not silently assume that the external accelerometer-
bias hold is eventually released.  Under MAG-CALL-SCHEDULE-v1 the internal
``accel_bias_locked_`` flag is forced clear within 10 s of gauged Live.  From
there the actual shipping language has two admissible continuations:

* no external hold: the same event enables A21 immediately;
* external hold: remain H18 for arbitrarily long; the first later hold release
  while Live applies the exact H->A variance-floor map.

An A21 execution may also be forced back to H18 by asserting the hold; the
literal shipping map zeros BA cross-covariances.  Thus eventual A21 is NOT a
universal premise of ALT.  Both held and active continuations stay in the
finite hybrid graph.

This module composes already proved control/covariance relations.  It does not
qualify the magnetic packet arithmetic, deployment clock roundoff, or storage.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_h18_a21_reachability as REACH


@dataclass(frozen=True)
class UnlockResult:
    state: CORE.State
    control: GATE.State
    live_to_unlock_upper: object
    mode_after_unlock: str
    held_continuation_allowed: bool


def _post_unlock_control(*, hold: bool, unlock_count: int = 250):
    # The reachability lemma proves that by this point first_time is finite and
    # the strict >1 s guard has fired.  Keep only the shipping control state
    # needed by subsequent setAccBiasHold events.
    return GATE.State(updates=unlock_count, first_time=F(1,25), locked=False, hold=hold)


def after_qualified_unlock(filter_state: CORE.State, cfg: GATE.Config, *, external_hold: bool):
    """Apply the mode consequence of the source-qualified unlock theorem.

    ``filter_state`` is the H18 MEKF state at the unlock boundary after the
    separately composed magnetic update.  No measurement acceptance premise is
    needed because the wrapper counter is attempt-based.
    """
    if not isinstance(filter_state, CORE.State) or filter_state.mode != 'H':
        raise ValueError('qualified unlock composition starts from H18 state')
    if not isinstance(cfg, GATE.Config):
        raise TypeError('shipping magnetic-bias config required')
    r = REACH.reach(external_hold=external_hold)
    control = _post_unlock_control(hold=external_hold, unlock_count=cfg.unlock_count)
    if external_hold:
        # Shipping clears the lock but does not enable BA while hold is active.
        return UnlockResult(filter_state, control, r.live_to_unlock_upper, 'H', True)
    # With no hold, the unlock event applies the same literal H->A floor as
    # set_acc_bias_updates_enabled(true).  Reuse set_hold by representing the
    # immediately-post-unlock control as held then releasing at the same Live
    # boundary; this calls the exact public finite edge and does not duplicate
    # its covariance algebra.
    pre = GATE.State(control.updates, control.first_time, False, True)
    edge = GATE.set_hold(pre, filter_state, cfg, hold=False, live=True)
    if edge.filter_state.mode != 'A' or not edge.enabled_bias_now:
        raise AssertionError('qualified no-hold unlock did not enter A21')
    return UnlockResult(edge.filter_state, control, r.live_to_unlock_upper, 'A', True)


def release_external_hold(result: UnlockResult, cfg: GATE.Config):
    """First hold release after lock-clear: H18 -> A21, or identity if active."""
    if not isinstance(result, UnlockResult) or not isinstance(cfg, GATE.Config):
        raise TypeError('unlock result and shipping config required')
    if result.state.mode == 'A':
        return result
    if not result.control.hold or result.control.locked:
        raise ValueError('release edge requires held H18 with internal lock clear')
    edge = GATE.set_hold(result.control, result.state, cfg, hold=False, live=True)
    if edge.filter_state.mode != 'A' or not edge.enabled_bias_now:
        raise AssertionError('hold release failed to enable A21')
    return UnlockResult(edge.filter_state, edge.state, result.live_to_unlock_upper, 'A', True)


def assert_external_hold(result: UnlockResult, cfg: GATE.Config):
    """Literal A21 -> H18 external-hold edge, retaining BA marginal."""
    if not isinstance(result, UnlockResult) or not isinstance(cfg, GATE.Config):
        raise TypeError('unlock result and shipping config required')
    edge = GATE.set_hold(result.control, result.state, cfg, hold=True, live=True)
    return UnlockResult(edge.filter_state, edge.state, result.live_to_unlock_upper,
                        edge.filter_state.mode, True)


def readiness():
    return {
      'qualified_schedule_to_internal_lock_clear_composed': True,
      'no_hold_unlock_to_A21_covariance_edge_composed': True,
      'indefinite_external_hold_H18_branch_retained': True,
      'post_unlock_hold_release_to_A21_edge_composed': True,
      'A21_to_H18_external_hold_edge_composed': True,
      'eventual_A21_not_assumed_for_arbitrary_hold_history': True,
      'both_H18_and_A21_continuations_in_hybrid_language': True,
      'deployment_clock_roundoff_closed': False,
      'complete_same_history_Live_word': False,
      'ALT_LIVE_PASS': False,
    }
