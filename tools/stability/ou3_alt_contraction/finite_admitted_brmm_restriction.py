"""Universal COMPLETE-BRMM history -> finite ALT restriction bridge.

The primary theorem quantifies over an *already admitted physical history*.
Admission is not something a 601-sample runtime object can infer from labels or
finite endpoint checks. This module makes that quantifier explicit and applies
``finite_complete_brmm_restriction`` to the fresh origin and every 5 ms step.

``AdmittedHistory`` is a theorem variable/hypothesis, not a runtime membership
certificate. No boolean "admitted" flag and no generator-token test is accepted.
``RestrictedOrigin`` and ``RestrictedSegment`` are theorem data obtained by
restricting that same quantified history. Necessary finite constraints are then
rechecked fail-closed by ``finite_source_continuation``.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_complete_brmm_restriction as RST
from tools.stability.ou3_alt_contraction import finite_source_continuation as SRC
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import bias_families as BIAS

QUALIFICATION='OU3_ALT_ADMITTED_PRIMARY_HISTORY_RESTRICTION_V1'
CANONICAL_SOURCE='COMPLETE_BRMM_NORMAL_LIVE_WORD'


def _restriction_definition():
    d=RST.build(); failures=RST.validate(d)
    if failures: raise RuntimeError('COMPLETE-BRMM restriction prerequisite failed: '+repr(failures))
    if d['canonical_source'] != CANONICAL_SOURCE or not d['primary_physics_is_definition']:
        raise RuntimeError('primary COMPLETE-BRMM source definition changed')
    if d['arbitrary_runtime_tokens_prove_primary_history_membership'] is not False:
        raise RuntimeError('runtime token membership must remain forbidden')
    return d


@dataclass(frozen=True)
class AdmittedHistory:
    """One arbitrary history under the theorem's COMPLETE-BRMM quantifier."""
    history_id: str
    definition_qualification: str = QUALIFICATION
    canonical_source: str = CANONICAL_SOURCE
    def __post_init__(self):
        if not isinstance(self.history_id,str) or not self.history_id:
            raise ValueError('quantified physical history identity required')
        if self.definition_qualification != QUALIFICATION or self.canonical_source != CANONICAL_SOURCE:
            raise ValueError('wrong admitted-history theorem definition')
        _restriction_definition()


@dataclass(frozen=True)
class RestrictedOrigin:
    """The t_L sample of the same quantified admitted history."""
    history: AdmittedHistory
    endpoint: PHYS.PhysicalKinematics
    def __post_init__(self):
        if not isinstance(self.history,AdmittedHistory) or not isinstance(self.endpoint,PHYS.PhysicalKinematics):
            raise TypeError('admitted history and fresh physical endpoint required')
        if self.endpoint.time != self.endpoint.live_origin:
            raise ValueError('admitted source origin must be the one-time Live origin')
        if any(self.endpoint.centered_S):
            raise ValueError('admitted source origin must have zero centered S')
        # CORE.Reference carries the non-qualifying ancestry labels. Require the
        # one used by the theorem bridge to agree with the quantified history.
        if getattr(self.endpoint,'history_id',None) != self.history.history_id:
            raise ValueError('fresh endpoint detached from quantified admitted history')


@dataclass(frozen=True)
class RestrictedSegment:
    """Exact kth restriction value of the quantified history (theorem datum)."""
    history: AdmittedHistory
    ordinal: int
    segment: PHYS.PhysicalSegment
    def __post_init__(self):
        if not isinstance(self.history,AdmittedHistory) or not isinstance(self.segment,PHYS.PhysicalSegment):
            raise TypeError('admitted history and exact physical segment required')
        if not isinstance(self.ordinal,int) or isinstance(self.ordinal,bool) or not 1 <= self.ordinal <= SRC.TRANSITIONS:
            raise ValueError('restriction ordinal must be one of the 600 transitions')
        expected=self.segment.before.live_origin+(self.ordinal-1)*SRC.DT
        if self.segment.before.time != expected or self.segment.after.time != expected+SRC.DT:
            raise ValueError('restricted segment clock detached from theorem sampling grid')
        for endpoint in (self.segment.before,self.segment.after):
            if getattr(endpoint,'history_id',self.history.history_id) != self.history.history_id:
                raise ValueError('restricted segment detached from quantified admitted history')


def source_root(history:AdmittedHistory, *, live_origin, bias_family:str):
    """Derive the finite-source root from one quantified admitted history."""
    if not isinstance(history,AdmittedHistory): raise TypeError('AdmittedHistory required')
    c=next((x for x in BIAS.contracts() if x.name==bias_family),None)
    if c is None: raise ValueError('unknown BIAS family')
    return SRC.SourceRoot(history.history_id,'primary-physical-history-restriction',F(live_origin),
                          c.name,c.parameter_token)


def qualify_origin(root:SRC.SourceRoot, restricted:RestrictedOrigin):
    if not isinstance(root,SRC.SourceRoot) or not isinstance(restricted,RestrictedOrigin):
        raise TypeError('finite source root and admitted-history origin required')
    if root.history_id != restricted.history.history_id:
        raise ValueError('finite source root detached from quantified admitted history')
    _restriction_definition()
    return SRC.origin_endpoint(root,restricted.endpoint)


def qualify_step(root:SRC.SourceRoot, restricted:RestrictedSegment,
                 witness:SRC.StepWitness):
    """Apply the universal restriction lemma, then finite necessary checks."""
    if not isinstance(root,SRC.SourceRoot) or not isinstance(restricted,RestrictedSegment):
        raise TypeError('finite source root and admitted-history restriction required')
    if root.history_id != restricted.history.history_id:
        raise ValueError('finite source root detached from quantified admitted history')
    if witness.ordinal != restricted.ordinal:
        raise ValueError('source witness ordinal detached from history restriction ordinal')
    _restriction_definition()
    return SRC.QualifiedPhysicalSegment(root,witness,restricted.segment)


def readiness():
    d=_restriction_definition()
    return {
      'universal_theorem_quantifier_over_admitted_primary_history_explicit':True,
      'finite_word_derived_by_restriction_not_runtime_membership_inference':True,
      'arbitrary_runtime_tokens_do_not_prove_COMPLETE_BRMM_membership':not d['arbitrary_runtime_tokens_prove_primary_history_membership'],
      'primary_history_restriction_implies_finite_physical_constraints':d['same_history_restriction_required'],
      'same_admitted_history_has_explicit_tL_origin_restriction':True,
      'one_same_history_identity_required_at_every_restricted_endpoint':True,
      'physical_state_coordinates_remain_free_of_proof_metadata':True,
      'restriction_ordinal_bound_to_canonical_5ms_grid':True,
      'BIAS_generating_history_attached_by_this_bridge':False,
      'sensor_disturbance_history_attached_by_this_bridge':False,
      'startup_sample_zero_equal_to_admitted_history_restriction_proved':False,
      'complete_600_step_shipping_word_composed_from_restrictions':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
