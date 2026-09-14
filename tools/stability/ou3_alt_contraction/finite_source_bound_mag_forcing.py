"""Expose the actual same-event magnetic disturbance supply for ALT.

``finite_live_magnetic_word`` already derives the MEKF magnetic forcing from the
same qualified physical field, current physical attitude, learned active
reference, applied hard-iron correction and raw magnetic residual.  This wrapper
makes those correlated terms explicit for the future ISS inequality instead of
replacing the effective residual by an independent disturbance box.

No contraction, smallness or deployment-roundoff claim is made here.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


@dataclass(frozen=True)
class MagForcing:
    rotated_reference_mismatch: tuple
    hard_iron_body: tuple
    applied_bias_correction: tuple
    raw_residual_body: tuple
    effective_residual: tuple

    def __post_init__(self):
        for name in ('rotated_reference_mismatch','hard_iron_body',
                     'applied_bias_correction','raw_residual_body','effective_residual'):
            object.__setattr__(self,name,tuple(M.vec(getattr(self,name),3)))
        expected=tuple(self.rotated_reference_mismatch[i]+self.hard_iron_body[i]
                       +self.applied_bias_correction[i]+self.raw_residual_body[i]
                       for i in range(3))
        if self.effective_residual != expected:
            raise ValueError('magnetic effective residual detached from same-event correlated components')

    @property
    def supply_vector(self):
        return (self.raw_residual_body + self.rotated_reference_mismatch
                + self.hard_iron_body + self.applied_bias_correction)


@dataclass(frozen=True)
class Result:
    word: object
    forcing: MagForcing | None

    @property
    def state(self): return self.word.state


def mag_step(state:WORD.State, **kwargs):
    """Execute the source-owned magnetic edge and retain its actual forcing."""
    out=WORD.mag_step(state,**kwargs)
    event=out.event.event
    # Outer-gated calls have no source packet and no MEKF magnetic forcing.
    if event.qualification is None:
        if event.effective_residual is not None:
            raise AssertionError('gated magnetic call unexpectedly produced effective residual')
        return Result(out,None)
    source=event.qualification.sample
    active=event.state.active
    memory=event.state.memory
    ref=event.filter.reference
    field_difference=tuple(source.model.world_field[i]-active.model.world_reference[i]
                           for i in range(3))
    rotated=tuple(SENSOR.q_rotate(ref.q_world_to_body,field_difference))
    correction=tuple(-x for x in memory.applied.total_bias)
    forcing=MagForcing(rotated,source.model.hard_iron_body,correction,
                       source.residual_body,event.effective_residual)
    return Result(out,forcing)


def readiness():
    return {
      'magnetic_effective_residual_derived_from_same_qualified_event':True,
      'raw_magnetic_residual_retained_separately':True,
      'active_reference_mismatch_retained_separately':True,
      'hard_iron_and_applied_correction_retained_separately':True,
      'independent_free_magnetic_innovation_forcing_forbidden_at_this_entry':True,
      'magnetic_forcing_smallness_claimed':False,
      'deployment_roundoff_supply_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
