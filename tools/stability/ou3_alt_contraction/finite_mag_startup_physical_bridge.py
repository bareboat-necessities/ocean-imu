"""Bridge startup magnetic source packets to the main finite physical endpoint.

``finite_mag_startup_source`` originally carried a small magnetic-only physical
endpoint.  The complete ALT word, however, already owns a richer
``PhysicalKinematics`` endpoint.  This bridge derives the magnetic endpoint from
that exact object so time and true attitude cannot diverge between the physical
word and startup north-learning path.

A separate history token is retained because ``PhysicalKinematics`` intentionally
contains kinematics only and does not itself certify COMPLETE-BRMM admission.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as MAG


def endpoint(physical:PHYS.PhysicalKinematics,history_id:str):
    if not isinstance(physical,PHYS.PhysicalKinematics):
        raise TypeError('main finite PhysicalKinematics endpoint required')
    if not isinstance(history_id,str) or not history_id:
        raise ValueError('persistent physical history token required')
    return MAG.PhysicalEndpoint(physical.time,physical.q_world_to_body,history_id)


def sample(physical:PHYS.PhysicalKinematics,history_id:str,model:MAG.Model,residual_body,packet_id):
    """Construct one startup magnetic source from the same finite physical endpoint."""
    return MAG.make_sample(endpoint(physical,history_id),model,residual_body,packet_id)


def assert_same_endpoint(physical:PHYS.PhysicalKinematics,mag_sample:MAG.Sample,history_id:str):
    if not isinstance(physical,PHYS.PhysicalKinematics) or not isinstance(mag_sample,MAG.Sample):
        raise TypeError('main finite endpoint and startup magnetic sample required')
    e=endpoint(physical,history_id)
    if mag_sample.physical!=e:
        raise ValueError('startup magnetic source detached from same finite physical endpoint')
    return True


def readiness():
    return {
      'startup_mag_time_from_main_PhysicalKinematics':True,
      'startup_mag_true_attitude_from_main_PhysicalKinematics':True,
      'duplicate_free_physical_endpoint_bridge':True,
      'history_token_retained_separately_from_kinematic_object':True,
      'COMPLETE_BRMM_admission_proved_by_object_identity':False,
      'magnetic_source_envelope_declared':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
