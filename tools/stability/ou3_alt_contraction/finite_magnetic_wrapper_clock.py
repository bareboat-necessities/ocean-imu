"""Exact outer-wrapper clock operands for ALT magnetic events.

Shipping has two distinct time coordinates at a magnetometer call:

* the physical/inner-MEKF time, advanced by the canonical IMU source; and
* ``SeaStateFusion_OU_III::t_``, an outer binary32 accumulator.

The existing finite magnetic algebra historically used the physical endpoint
clock for both.  This module materializes the correct outer clock on the
canonical 5 ms grid, including the literal shipping fallback rule used by
continuous hard-iron accumulation/application:

    dt = current-last if current > last else mag_sample_dt_sec.

It does not silently promote the existing magnetic composer: until that composer
consumes these operands for its outer delay/refinement/continuous branches,
``dual_clock_magnetic_word_composed`` remains false.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

QUALIFICATION='OU3_ALT_MAGNETIC_WRAPPER_CLOCK_V1'
DEFAULT_SAMPLE_DT=F(1,200)
DEFAULT_MAG_DELAY=F(7)
DEFAULT_REFINEMENT_START=F(90)


def R(x): return M.rational(x)


@dataclass(frozen=True)
class Timestamp:
    physical_time:F
    source_step:int
    wrapper_time:F
    def __post_init__(self):
        p=R(self.physical_time); w=R(self.wrapper_time)
        if not isinstance(self.source_step,int) or isinstance(self.source_step,bool):
            raise TypeError('canonical source step integer required')
        if CLOCK.canonical_step_index(p)!=self.source_step:
            raise ValueError('physical time detached from canonical 5 ms source step')
        expected=CLOCK.clock_at_step(self.source_step)
        if w!=expected:
            raise ValueError('wrapper timestamp detached from exact binary32 recurrence')
        object.__setattr__(self,'physical_time',p); object.__setattr__(self,'wrapper_time',w)


def at_physical_time(physical_time):
    p=R(physical_time); k=CLOCK.canonical_step_index(p)
    return Timestamp(p,k,CLOCK.clock_at_step(k))


def shipping_elapsed(current:Timestamp, previous_wrapper_time, *, fallback_dt=DEFAULT_SAMPLE_DT):
    """Literal outer-wrapper dt selection used by hard-iron sample/apply clocks."""
    if not isinstance(current,Timestamp): raise TypeError('certified wrapper timestamp required')
    fallback=R(fallback_dt)
    if fallback<=0: raise ValueError('positive configured magnetic fallback dt required')
    if previous_wrapper_time is None: return fallback
    prev=R(previous_wrapper_time)
    return current.wrapper_time-prev if current.wrapper_time>prev else fallback


def outer_mag_delay_passed(ts:Timestamp, *, delay=DEFAULT_MAG_DELAY):
    return ts.wrapper_time>=R(delay)


def refinement_due(ts:Timestamp, *, start=DEFAULT_REFINEMENT_START):
    return ts.wrapper_time>=R(start)


def readiness():
    c=CLOCK.readiness()
    return {
      'qualification':QUALIFICATION,
      'physical_and_outer_wrapper_clocks_are_distinct_coordinates':True,
      'exact_binary32_wrapper_timestamp_derived_from_canonical_physical_endpoint':c['exact_wrapper_clock_lookup_available_on_certified_grid'],
      'shipping_nonadvancing_clock_fallback_dt_materialized':True,
      'outer_mag_delay_predicate_uses_wrapper_clock':True,
      'outer_refinement_start_predicate_uses_wrapper_clock':True,
      'continuous_sample_and_apply_dt_use_wrapper_clock':True,
      'canonical_prefix_wrapper_clock_arithmetic_closed':c['canonical_5ms_wrapper_clock_prefix_binary32_closed'],
      'dual_clock_magnetic_word_composed':False,
      'indefinite_wrapper_clock_lifetime_closed':False,
      'deployment_exp_solver_roundoff_closed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
