"""Declared deployment scope for the independent OU-III ALT theorem.

The shipping filter supports an optional steady-wind-heel retargeting feature.
ALT does not certify that optional branch.  The certified deployment language is
restricted to histories for which

    wind_heel_rad_ == 0

from construction onward and ``update_wind_heel`` is never invoked.  This is a
scope restriction, not a shipping-code change: the shipping member defaults to
zero.  Consequently the virtual B' frame equals the physical body frame B,
``deheel_vector_`` is the identity, ``quaternion_boat``/``set_quaternion_boat``
contain no nontrivial heel rotation, and no body-frame retarget event belongs to
the certified hybrid language.

A future proof may widen this scope by proving the retarget map.  Until then an
execution with nonzero heel or any wind-heel update is outside ALT.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

QUALIFICATION='OU3_ALT_DEPLOYMENT_SCOPE_ZERO_WIND_HEEL_V1'

@dataclass(frozen=True)
class Scope:
    wind_heel_rad:F=F(0)
    update_wind_heel_calls:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        h=F(self.wind_heel_rad)
        if h != 0:
            raise ValueError('ALT certified deployment excludes nonzero wind heel')
        if int(self.update_wind_heel_calls) != 0:
            raise ValueError('ALT certified deployment excludes update_wind_heel calls')
        if self.qualification != QUALIFICATION:
            raise ValueError('wrong ALT deployment-scope qualification')
        object.__setattr__(self,'wind_heel_rad',h)
        object.__setattr__(self,'update_wind_heel_calls',0)


def certified_scope():
    return Scope()


def assert_certified_scope(scope:Scope):
    if not isinstance(scope,Scope):
        raise TypeError('ALT zero-wind-heel deployment scope required')
    if scope.wind_heel_rad != 0 or scope.update_wind_heel_calls != 0:
        raise ValueError('wind-heel branch is outside ALT theorem scope')
    return True


def readiness():
    return {
      'shipping_default_wind_heel_zero':True,
      'wind_heel_feature_excluded_from_ALT':True,
      'update_wind_heel_event_excluded_from_hybrid_language':True,
      'body_prime_equals_physical_body_in_certified_scope':True,
      'shipping_filter_changed':False,
      'original_proof_track_changed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
