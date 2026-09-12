"""Deterministic magnetic source qualification for the ALT theorem.

The finite startup magnetic graph retains the exact source identity

    m_raw_B = R_true B_W + b_HI + n_m.

MAG-BMM150-DET-v1 is an engineering admission class for a commissioned
BMM150 installation, not a Bosch guarantee for arbitrary installations:

    20 uT <= ||B_W||_2 <= 75 uT,
    ||(B_Wx,B_Wy)||_2 >= 15 uT,
    ||b_HI||_2 <= 5 uT,
    ||n_m||_2 <= 2 uT per accepted sample.

The Earth-field interval encloses the NOAA/WMM surface total-field range with
margin.  The hard-iron and residual limits are commissioned-installation/source
requirements.  They are intentionally deterministic and materially wider than
BMM150 RMS output noise, but tight enough to preserve a useful horizontal-north
capture margin under the already-declared 0.02 rad startup tilt-direction error.
Samples outside them are outside the theorem even when the sensor itself is far
from electrical saturation.  The horizontal lower bound is required for
deterministic yaw observability; a total-field bound alone cannot prove north
capture near a magnetic pole.

Rmag remains a stochastic/model covariance and is never used as a deterministic
source bound.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as SOURCE


def R(x): return M.rational(x)

ASSUMPTION_ID='MAG-BMM150-DET-v1'
DEFAULT_WORLD_FIELD_NORM_MIN=F(20)
DEFAULT_WORLD_FIELD_NORM_MAX=F(75)
DEFAULT_WORLD_FIELD_HORIZONTAL_MIN=F(15)
DEFAULT_HARD_IRON_NORM_MAX=F(5)
DEFAULT_RESIDUAL_NORM_MAX=F(2)


@dataclass(frozen=True)
class Envelope:
    world_field_norm_min:F=DEFAULT_WORLD_FIELD_NORM_MIN
    world_field_norm_max:F=DEFAULT_WORLD_FIELD_NORM_MAX
    world_field_horizontal_min:F=DEFAULT_WORLD_FIELD_HORIZONTAL_MIN
    hard_iron_norm_max:F=DEFAULT_HARD_IRON_NORM_MAX
    residual_norm_max:F=DEFAULT_RESIDUAL_NORM_MAX
    assumption_id:str=ASSUMPTION_ID
    def __post_init__(self):
        for n in ('world_field_norm_min','world_field_norm_max','world_field_horizontal_min',
                  'hard_iron_norm_max','residual_norm_max'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative deterministic magnetic source envelope required')
            object.__setattr__(self,n,v)
        if self.world_field_norm_min>self.world_field_norm_max:
            raise ValueError('magnetic field lower bound cannot exceed upper bound')
        if self.world_field_horizontal_min>self.world_field_norm_max:
            raise ValueError('horizontal field lower bound cannot exceed total-field upper bound')
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
class HorizontalNormWitness:
    xy:tuple
    norm:F
    def __post_init__(self):
        if len(self.xy)!=2: raise ValueError('two horizontal field components required')
        xy=tuple(R(x) for x in self.xy); n=R(self.norm)
        if n<0 or n*n!=xy[0]*xy[0]+xy[1]*xy[1]:
            raise ValueError('exact horizontal-field norm witness required')
        object.__setattr__(self,'xy',xy); object.__setattr__(self,'norm',n)


@dataclass(frozen=True)
class Qualification:
    sample:SOURCE.Sample
    envelope:Envelope
    qualified:bool


def default_envelope():
    return Envelope()


def qualify(sample:SOURCE.Sample,envelope:Envelope|None=None,*,field_norm:NormWitness|None=None,
            horizontal_norm:HorizontalNormWitness|None=None,
            hard_iron_norm:NormWitness|None=None,residual_norm:NormWitness|None=None):
    if not isinstance(sample,SOURCE.Sample): raise TypeError('startup magnetic source sample required')
    envelope=default_envelope() if envelope is None else envelope
    if not isinstance(envelope,Envelope): raise TypeError('magnetic source envelope required')
    if envelope.assumption_id!=ASSUMPTION_ID:
        raise ValueError('ALT magnetic source promotion requires MAG-BMM150-DET-v1')

    if not isinstance(field_norm,NormWitness) or field_norm.vector!=tuple(sample.model.world_field):
        raise ValueError('world field norm witness detached from same source sample')
    if not (envelope.world_field_norm_min<=field_norm.norm<=envelope.world_field_norm_max):
        raise ValueError('world field violates declared deterministic source envelope')

    xy=tuple(sample.model.world_field[:2])
    if not isinstance(horizontal_norm,HorizontalNormWitness) or horizontal_norm.xy!=xy:
        raise ValueError('horizontal field norm witness detached from same source sample')
    if horizontal_norm.norm<envelope.world_field_horizontal_min:
        raise ValueError('horizontal field violates deterministic yaw-observability envelope')

    triples=(
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
    if q.envelope.assumption_id!=ASSUMPTION_ID:
        raise ValueError('wrong magnetic theorem-source assumption')
    return True


def readiness():
    return {
      'magnetic_source_identity_is_separate_from_source_envelope':True,
      'Rmag_is_not_used_as_deterministic_noise_bound':True,
      'named_BMM150_deterministic_envelope_declared':True,
      'horizontal_yaw_observability_lower_bound_declared':True,
      'canonical_theorem_currently_declares_magnetic_envelope':True,
      'individual_magnetic_sample_can_be_source_qualified':True,
      'source_uniform_magnetic_word_qualified':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
