"""Quantified admitted BIAS0/1/2 history for the source-owned ALT word.

``finite_bias_history_restriction`` proves the family-level restriction theorem.
This module supplies the corresponding theorem variable: one already-admitted
physical bias history, one fixed recurrence factor over the 5 ms word, and the
kth restriction datum of that same history.  It is not a runtime membership
oracle and does not infer hardware admission.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_bias_history_restriction as BHR
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_BIAS_HISTORY_V1'


def _family(name):
    d=BHR.build(); failures=BHR.validate(d)
    if failures: raise RuntimeError('BIAS restriction theorem invalid: '+repr(failures))
    for r in BHR.restrictions():
        if r.name==name: return r
    raise ValueError('unknown BIAS family')


@dataclass(frozen=True)
class AdmittedBiasHistory:
    history_id: str
    family: str
    phi_true: F
    qualification: str = QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.history_id,str) or not self.history_id:
            raise ValueError('quantified bias-history identity required')
        if self.qualification != QUALIFICATION:
            raise ValueError('wrong admitted-bias theorem definition')
        r=_family(self.family); phi=M.rational(self.phi_true)
        lo=F.from_float(r.phi_lo); hi=F.from_float(r.phi_hi)
        if not lo <= phi <= hi:
            raise ValueError('physical bias recurrence factor outside admitted family')
        if self.family=='BIAS2' and phi != 1:
            raise ValueError('canonical BIAS2 restriction uses non-relaxing phi=1')
        object.__setattr__(self,'phi_true',phi)


@dataclass(frozen=True)
class RestrictedBiasStep:
    history: AdmittedBiasHistory
    ordinal: int
    segment: PHYS.PhysicalSegment
    def __post_init__(self):
        if not isinstance(self.history,AdmittedBiasHistory) or not isinstance(self.segment,PHYS.PhysicalSegment):
            raise TypeError('admitted bias history and physical segment required')
        if not isinstance(self.ordinal,int) or isinstance(self.ordinal,bool) or not 1 <= self.ordinal <= SOURCE.TRANSITIONS:
            raise ValueError('bias restriction ordinal must be one of 600 transitions')
        if M.rational(self.segment.phi_true) != self.history.phi_true:
            raise ValueError('segment bias factor detached from one admitted bias history')
        expected=tuple(self.history.phi_true*self.segment.before.beta[i]+self.segment.bias_driver[i]
                       for i in range(3))
        if tuple(self.segment.after.beta) != expected:
            raise ValueError('bias restriction violates same-history beta recurrence')
        # Family hard bounds are checked again by QualifiedPhysicalSegment at
        # composition; this object establishes ancestry/recurrence, not a
        # second independent source envelope.


def readiness():
    d=BHR.build(); failures=BHR.validate(d)
    if failures: raise RuntimeError('BIAS restriction theorem invalid: '+repr(failures))
    return {
      'BIAS0_BIAS1_BIAS2_admitted_history_quantifier_available':True,
      'one_fixed_physical_phi_carried_per_BIAS_history':True,
      'BIAS2_nonrelaxing_phi_one_retained':d['BIAS2_canonical_non_relaxing_phi_is_one'],
      'same_history_beta_driver_recurrence_checked_each_restriction':True,
      'runtime_tokens_do_not_establish_bias_admission':True,
      'assembled_sensor_hardware_admission_inferred':False,
      'binary32_bias_runtime_correspondence_closed':False,
      'ALT_LIVE_PASS':False,
    }
