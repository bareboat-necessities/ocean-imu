"""Bind the fresh startup Live truth to sample zero of an admitted history.

The theorem quantifies over an already admitted COMPLETE-BRMM physical history.
A finite endpoint cannot prove that admission, but the theorem still needs an
explicit restriction datum saying which value that admitted history has at the
one-time Live origin.  This module carries exactly that datum and requires the
fresh joint24 Reference produced by startup composition to be that same object.

This closes an identity/ancestry edge only.  It does not infer admission from
endpoint checks, prove startup reachability for every history, or enclose any
deployment floating-point arithmetic.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ADMIT
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_brmm_moment_prefix as MOMENTS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as LIVE


@dataclass(frozen=True)
class RestrictedOrigin:
    """Exact sample-zero value of one quantified admitted physical history."""
    history: ADMIT.AdmittedHistory
    reference: CORE.Reference
    def __post_init__(self):
        if not isinstance(self.history,ADMIT.AdmittedHistory):
            raise TypeError('admitted COMPLETE-BRMM history required')
        if not isinstance(self.reference,CORE.Reference):
            raise TypeError('fresh physical Reference required')
        r=self.reference
        if r.time != r.live_origin:
            raise ValueError('admitted sample zero must occur at the one-time Live origin')
        if any(r.centered_S):
            raise ValueError('admitted sample zero must have centered S = 0')
        # These are necessary consequences of the primary physical definition,
        # not a membership test for the history quantified above.
        MOMENTS.check_endpoint(r)


@dataclass(frozen=True)
class State:
    admitted: LIVE.State
    origin: RestrictedOrigin
    def __post_init__(self):
        if not isinstance(self.admitted,LIVE.State) or not isinstance(self.origin,RestrictedOrigin):
            raise TypeError('admitted Live state and admitted sample-zero restriction required')
        if self.origin.history != self.admitted.admitted_history:
            raise ValueError('sample-zero restriction detached from carried admitted history')
        ref=self.admitted.live_word.live.live.live.mekf.reference
        if ref != self.origin.reference:
            raise ValueError('fresh Live physical Reference is not admitted-history sample zero')
        root=self.admitted.live_word.source.root
        if root.history_id != self.origin.history.history_id:
            raise ValueError('finite source root detached from admitted sample-zero history')
        if root.live_origin != self.origin.reference.live_origin:
            raise ValueError('finite source Live origin detached from admitted sample zero')
        if self.admitted.live_word.source.steps:
            raise ValueError('sample-zero binding applies before the first IMU transition')


def bind(admitted:LIVE.State, origin:RestrictedOrigin):
    """Attach theorem sample-zero restriction to the exact fresh Live truth."""
    return State(admitted,origin)


def readiness():
    return {
      'admitted_history_sample_zero_restriction_datum_explicit':True,
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
