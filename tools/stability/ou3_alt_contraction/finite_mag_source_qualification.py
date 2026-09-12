"""Fail-closed magnetic source qualification guard for the ALT theorem.

The finite startup magnetic graph now retains the exact source identity

    m_raw_B = R_true B_W + b_HI + n_m.

That is not yet a deterministic theorem source class.  The current canonical
COMPLETE-BRMM/BIAS assumptions do not declare finite envelopes for the magnetic
world field, body-fixed hard iron, or magnetometer residual.  This module makes
that omission explicit: local finite identities may use the source graph, but a
source-qualified finite master cannot consume it until named deterministic
bounds are supplied by theorem/source assumptions.

No numerical values are inferred from simulation, sensor datasheets, Rmag, or
sample statistics here.  Measurement covariance is not a deterministic source
bound.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as SOURCE


def R(x): return M.rational(x)

@dataclass(frozen=True)
class Envelope:
    world_field_norm_max:F
    hard_iron_norm_max:F
    residual_norm_max:F
    assumption_id:str
    def __post_init__(self):
        for n in ('world_field_norm_max','hard_iron_norm_max','residual_norm_max'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative deterministic magnetic source envelope required')
            object.__setattr__(self,n,v)
        if not isinstance(self.assumption_id,str) or not self.assumption_id:
            raise ValueError('named theorem/source assumption id required')

@dataclass(frozen=True)
class NormWitness:
    vector:tuple
    norm:F
    def __post_init__(self):
        v=tuple(M.vec(self.vector,3)); n=R(self.norm)
        if n<0 or n*n!=M.dot(v,v): raise ValueError('exact vector norm witness required')
        object.__setattr__(self,'vector',v); object.__setattr__(self,'norm',n)

@dataclass(frozen=True)
class Qualification:
    sample:SOURCE.Sample
    envelope:Envelope
    qualified:bool


def qualify(sample:SOURCE.Sample,envelope:Envelope|None,*,field_norm:NormWitness|None=None,
            hard_iron_norm:NormWitness|None=None,residual_norm:NormWitness|None=None):
    if not isinstance(sample,SOURCE.Sample): raise TypeError('startup magnetic source sample required')
    if envelope is None:
        raise ValueError('magnetic source envelopes are undeclared; theorem source qualification must fail closed')
    if not isinstance(envelope,Envelope): raise TypeError('magnetic source envelope required')
    triples=(
      (field_norm,sample.model.world_field,envelope.world_field_norm_max,'world field'),
      (hard_iron_norm,sample.model.hard_iron_body,envelope.hard_iron_norm_max,'hard iron'),
      (residual_norm,sample.residual_body,envelope.residual_norm_max,'mag residual'),
    )
    for witness,vec,bound,label in triples:
        if not isinstance(witness,NormWitness) or witness.vector!=tuple(vec):
            raise ValueError(f'{label} norm witness detached from same source sample')
        if witness.norm>bound:
            raise ValueError(f'{label} exceeds declared deterministic source envelope')
    return Qualification(sample,envelope,True)


def assert_source_qualified(q:Qualification|None):
    if not isinstance(q,Qualification) or not q.qualified:
        raise ValueError('unqualified magnetic source cannot enter finite master/storage boundary')
    return True


def readiness():
    return {
      'magnetic_source_identity_is_separate_from_source_envelope':True,
      'Rmag_is_not_used_as_deterministic_noise_bound':True,
      'missing_magnetic_envelope_fails_closed':True,
      'named_assumption_required_before_source_promotion':True,
      'canonical_theorem_currently_declares_magnetic_envelope':False,
      'source_uniform_magnetic_word_qualified':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
