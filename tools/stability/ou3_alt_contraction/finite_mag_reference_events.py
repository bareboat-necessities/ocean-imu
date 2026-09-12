"""Named magnetic-reference write events for the finite ALT word.

Shipping has one write funnel, setMagWorldRef_, used by the startup provisional
reference, second-stage refinement and continuous hard-iron correction. Ordinary
updateMag never changes the world reference. This module turns that fact into a
persistent generation ledger without claiming the estimator-generated vectors
are source-qualified yet.

A write witness is therefore not proof of the numeric reference; it names which
shipping producer supplied a value and forces generation/root continuity. The
producer-specific finite arithmetic remains open and fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_mag_reference_runtime as REF

REASONS=('startup_provisional','live_refinement','continuous_hard_iron')

@dataclass(frozen=True)
class WriteWitness:
    reason:str
    model:MAG.Model
    producer_root:str
    def __post_init__(self):
        if self.reason not in REASONS: raise ValueError('declared magnetic reference write reason required')
        if not isinstance(self.model,MAG.Model): raise TypeError('written magnetic model required')
        if not isinstance(self.producer_root,str) or not self.producer_root:
            raise ValueError('persistent magnetic producer root required')


def first_write(witness:WriteWitness):
    if not isinstance(witness,WriteWitness): raise TypeError('mag reference write witness required')
    if witness.reason!='startup_provisional':
        raise ValueError('first active magnetic reference must be startup provisional write')
    return REF.State(witness.model,0,witness.producer_root)


def rewrite(active:REF.State,witness:WriteWitness):
    if not isinstance(active,REF.State) or not isinstance(witness,WriteWitness):
        raise TypeError('active reference and write witness required')
    if witness.producer_root!=active.root:
        raise ValueError('magnetic reference producer root cannot restart')
    if witness.reason=='startup_provisional':
        raise ValueError('startup provisional reference cannot restart after active generation exists')
    # Rmag comes from constructor sigma_m and is not rewritten by setMagWorldRef_.
    if witness.model.sigma_internal!=active.model.sigma_internal:
        raise ValueError('reference write cannot change constructor Rmag/sigma_m')
    return REF.State(witness.model,active.generation+1,active.root)


def no_write(active:REF.State):
    if not isinstance(active,REF.State): raise TypeError('active magnetic reference required')
    return active


def readiness():
    return {
      'all_mag_reference_writes_use_named_setMagWorldRef_funnel':True,
      'startup_provisional_is_generation_zero':True,
      'refinement_and_continuous_writes_increment_generation':True,
      'reference_write_preserves_constructor_Rmag_sigma':True,
      'ordinary_updateMag_preserves_generation':True,
      'MagAutoTuner_startup_numeric_output_attached':False,
      'MagAutoTuner_refinement_numeric_output_attached':False,
      'continuous_hard_iron_numeric_output_attached':False,
      'mag_reference_validity_binary32_gates_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
