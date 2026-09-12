"""Persistent active magnetic model identity for the finite ALT word.

The shipping wrapper stores the world magnetic reference that it writes into the
MEKF; updateMag does not accept a new world reference.  Therefore each async
measurement must inherit one persistent active model rather than choosing a free
B_world/Rmag pair per event.

This layer only enforces persistence.  How startup MagAutoTuner, refinement and
continuous hard-iron correction create each successive active reference remains
an explicit source/frontend obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS

@dataclass(frozen=True)
class State:
    model:MAG.Model
    generation:int
    root:str
    def __post_init__(self):
        if not isinstance(self.model,MAG.Model): raise TypeError('active magnetic model required')
        if not isinstance(self.generation,int) or self.generation<0: raise ValueError('nonnegative magnetic reference generation required')
        if not isinstance(self.root,str) or not self.root: raise ValueError('persistent magnetic reference root required')


def sample(active:State,physical:PHYS.PhysicalKinematics,raw_body,residual_internal,*,deheel=None):
    if not isinstance(active,State) or not isinstance(physical,PHYS.PhysicalKinematics):
        raise TypeError('active magnetic reference and physical endpoint required')
    kw={} if deheel is None else {'deheel_body_to_internal':deheel}
    return MAG.Sample(physical,raw_body,residual_internal,active.model,**kw)


def require_same(active:State,packet:MAG.Sample):
    if not isinstance(active,State) or not isinstance(packet,MAG.Sample):
        raise TypeError('active reference and magnetic packet required')
    if packet.model != active.model:
        raise ValueError('updateMag packet detached from persistent active magnetic model')
    return packet


def readiness():
    return {
      'one_persistent_world_reference_and_Rmag_per_active_generation':True,
      'updateMag_cannot_choose_free_world_reference':True,
      'updateMag_cannot_choose_free_Rmag':True,
      'startup_mag_reference_generation_attached':False,
      'second_stage_mag_refinement_attached':False,
      'continuous_hard_iron_reference_update_attached':False,
      'configured_sigma_m_constructor_ancestry_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
