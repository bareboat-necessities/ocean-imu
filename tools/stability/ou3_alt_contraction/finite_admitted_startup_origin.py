"""Bind fresh startup Live truth to the canonical admitted-history sample zero.

``finite_admitted_brmm_restriction.RestrictedOrigin`` is the unique theorem
datum for the value of an already-admitted COMPLETE-BRMM history at t_L.  This
module connects that datum to the exact fresh joint24 physical Reference.  It
never infers history admission from endpoint checks.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ADMIT
from tools.stability.ou3_alt_contraction import finite_brmm_moment_prefix as MOMENTS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as LIVE

RestrictedOrigin = ADMIT.RestrictedOrigin


@dataclass(frozen=True)
class State:
    admitted: LIVE.State
    origin: RestrictedOrigin
    def __post_init__(self):
        if not isinstance(self.admitted,LIVE.State) or not isinstance(self.origin,RestrictedOrigin):
            raise TypeError('admitted Live state and canonical admitted origin required')
        if self.origin.history != self.admitted.admitted_history:
            raise ValueError('sample-zero restriction detached from carried admitted history')
        ref=self.admitted.live_word.live.live.live.mekf.reference
        if ref != self.origin.endpoint:
            raise ValueError('fresh Live physical Reference is not admitted-history sample zero')
        root=self.admitted.live_word.source.root
        qualified=ADMIT.qualify_origin(root,self.origin)
        if qualified.endpoint != ref:
            raise AssertionError('canonical admitted origin qualification changed endpoint')
        if root.live_origin != self.origin.endpoint.live_origin:
            raise ValueError('finite source Live origin detached from admitted sample zero')
        if self.admitted.live_word.source.steps:
            raise ValueError('sample-zero binding applies before the first IMU transition')
        MOMENTS.check_endpoint(ref)


def bind(admitted:LIVE.State, origin:RestrictedOrigin):
    return State(admitted,origin)


def readiness():
    a=ADMIT.readiness()
    return {
      'canonical_admitted_history_sample_zero_restriction_consumed':a['same_admitted_history_has_explicit_tL_origin_restriction'],
      'fresh_joint24_physical_reference_equals_admitted_sample_zero':True,
      'same_one_time_Live_origin_shared_by_source_root_and_admitted_origin':True,
      'centered_S_zero_shared_at_fresh_admitted_origin':True,
      'finite_endpoint_checks_used_as_membership_oracle':False,
      'startup_sample_zero_equal_to_admitted_history_restriction_proved':True,
      'startup_reachability_for_every_admitted_history_proved':False,
      'deployment_binary32_correspondence_closed':False,
      'complete_same_history_startup_to_600_step_shipping_word':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
