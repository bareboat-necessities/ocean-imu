"""Fail-closed deployment obstruction for an *indefinite* ALT theorem.

The mathematical contraction target is indefinite, but ALT also targets the
current shipping implementation and finite-precision execution. Two independent
runtime-lifetime obligations prevent an all-time machine-execution theorem from
being claimed at present:

1. the outer wrapper clock is binary32.  The exact recurrence module exhibits
   RN32(2^17 + RN32(0.005)) == 2^17, so the canonical 5 ms increment eventually
   ceases to advance that shipping clock;
2. ``mag_updates_applied_`` is a signed ``int`` incremented on every post-delay
   magnetometer call.  The current MAG-CALL-SCHEDULE-v1 gives only upper bounds
   on call gaps, hence no uniform call-count upper bound and no signed-overflow
   proof, even on a fixed compact window.

This is NOT a dynamical-instability counterexample.  It is a deployment-language
obstruction: repeated finite words cannot be promoted to an indefinite theorem
until the shipping arithmetic is proved safe, the theorem declares a finite
lifetime within the proved arithmetic horizon, or the implementation is changed
under a separately authorized scope.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as MAG

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_INDEFINITE_DEPLOYMENT_OBSTRUCTION_V1'


@dataclass(frozen=True)
class Report:
    wrapper_clock_finite_prefix_closed: bool
    wrapper_clock_indefinite_closed: bool
    exact_clock_stall_witness_present: bool
    signed_mag_counter_increment_present: bool
    signed_mag_counter_saturation_present: bool
    mag_schedule_uniform_call_upper_present: bool
    signed_mag_counter_lifetime_closed: bool
    indefinite_current_shipping_execution_closed: bool


def _counter_source_shape():
    s=SOURCE.read_text()
    increment=s.count('mag_updates_applied_++;')
    declaration=s.count('int  mag_updates_applied_ = 0;')
    # Saturating/clamped increment would need to replace the literal unchecked
    # post-measurement increment before this proof can become true.
    saturation=('std::numeric_limits<int>::max()' in s and 'mag_updates_applied_' in s)
    return declaration>=1 and increment>=1, saturation


def build():
    c=CLOCK.readiness(); m=MAG.counter_lifetime(); inc,sat=_counter_source_shape()
    indefinite=bool(c['indefinite_wrapper_clock_lifetime_closed'] and
                    m.no_signed_overflow_proved and inc and sat)
    return Report(
        bool(c['canonical_5ms_wrapper_clock_prefix_binary32_closed']),
        bool(c['indefinite_wrapper_clock_lifetime_closed']),
        bool(c['late_time_2pow17_stall_witness_present']),
        inc,sat,
        m.uniform_call_count_upper is not None,
        bool(m.no_signed_overflow_proved),
        indefinite)


def assert_indefinite_shipping_execution_closed(report:Report|None=None):
    r=build() if report is None else report
    if not isinstance(r,Report) or not r.indefinite_current_shipping_execution_closed:
        raise RuntimeError('indefinite current-shipping execution is not deployment-arithmetic closed')
    return True


def readiness():
    r=build()
    return {
      'qualification':QUALIFICATION,
      'canonical_finite_wrapper_clock_prefix_closed':r.wrapper_clock_finite_prefix_closed,
      'exact_binary32_late_time_clock_stall_witness_present':r.exact_clock_stall_witness_present,
      'shipping_signed_mag_counter_unchecked_increment_present':r.signed_mag_counter_increment_present,
      'shipping_signed_mag_counter_saturation_present':r.signed_mag_counter_saturation_present,
      'mag_schedule_uniform_call_count_upper_present':r.mag_schedule_uniform_call_upper_present,
      'shipping_signed_mag_counter_lifetime_closed':r.signed_mag_counter_lifetime_closed,
      'indefinite_wrapper_clock_lifetime_closed':r.wrapper_clock_indefinite_closed,
      'indefinite_current_shipping_execution_closed':r.indefinite_current_shipping_execution_closed,
      'finite_word_contraction_invalidated_by_this_obstruction':False,
      'dynamical_instability_claimed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
