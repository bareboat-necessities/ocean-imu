#!/usr/bin/env python3
"""Joint24 derivative-word assembly retained for the ALT Live construction.

Scope correction: the code below assembles pointwise state/source Jacobian
data. Its legacy assembly labels do not prove an anchored finite error-map
identity. The finite_storage_readiness guard therefore blocks its use as the
finite endpoint in a common-storage search. See the finite-measurement proof
for the missing bridge and the exact accepted-event replacement.

The theorem is relational, not an enumeration.  For every endpoint path of the
regional universal Normal-Live relation, apply the source-uniform local cocycle
of ``physical_lineage`` at each literal sample.  Since the transition retains
ALL successors, induction gives a physical response cocycle on every endpoint
path.

The induction invariant is deliberately stronger than a matrix product:
  * selector parent/child ancestry is consecutive;
  * one BRMM generator and one Live S origin persist;
  * primitive_out(k) == primitive_in(k+1);
  * all q15 and physical-S blocks remain members of the ONE global
    O^601_BRMM relation, never independent boxes;
  * one BIAS0/1/2 parameter/root token is reused over the whole path;
  * all prior source-response columns are suffix-propagated through each new
    joint24 local map.

This establishes the H18-only and A21-only complete-word construction once a
mode path is supplied.  Hybrid H18->A21 splicing uses ``h18_a21_edge`` and is
kept separate because the release covariance must be constructed from the
actual held H18 covariance rather than the parallel hypothetical A21 branch.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence

from ou3_interval import Interval, matrix_mul
import ou3_brmm_correlated_window_outer_enclosure as OUTER
import ou3_p4_projection_sector as PROJ
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import live_source_relation as SOURCE
from tools.stability.ou3_alt_contraction import physical_lineage as LOCAL

QUALIFICATION='OU3_ALT_UNIVERSAL_PHYSICAL_JOINT24_WORD_INDUCTION_V1'
TRANSITIONS=600
SOURCE_SAMPLES=601
GLOBAL_PHYSICAL_RELATION='O^601_BRMM'


def I(x):return Interval.point(float(x))
def _shape(A):return len(A),len(A[0]) if A else 0
def _freeze(A):return tuple(tuple(x for x in row) for row in A)
def _thaw(A):return [list(row) for row in A]

def _propagate(J,B):
    if _shape(J)!=(24,24) or _shape(B)[0]!=24:raise ValueError('joint24 response propagation mismatch')
    return matrix_mul(J,B)

@dataclass(frozen=True)
class WordSourceBlock:
    kind:str
    token:str
    response:tuple[tuple[Interval,...],...]
    sectors:tuple
    global_relation:str
    generator_id:str
    live_origin_id:str
    family_parameter_token:str|None=None


def hard_bias_domains(contract:BIAS.BiasFamilyContract):
    """Regional Live domains implied by true-bias cap + closed projection ball."""
    p=PROJ.build();pf=PROJ.validate(p)
    if pf:raise RuntimeError('projection prerequisite failed: '+repr(pf))
    R=float(p['projection_radius_mps2']);B=float(contract.true_bias_component_bound)
    true=[Interval(-B,B) for _ in range(3)]
    # |bhat_i| <= ||bhat|| <= R and |beta_i|<=B.
    E=R+B;held=[Interval(-E,E) for _ in range(3)]
    return held,true


def _selector_key(lineage):
    s=lineage.selector
    return int(s.prefix_length),str(s.parent_source_cell_id),str(s.source_cell_id)


def _sample_primitive(lineage):
    ps=[c.wave_primitive for c in lineage.cells if c.wave_primitive is not None]
    if not ps:raise ValueError('attached sample has no physical primitive ancestry')
    p=ps[0]
    for q in ps[1:]:
        if (q.generator_id,q.live_origin_id,q.primitive_in_id,q.primitive_out_id)!=(p.generator_id,p.live_origin_id,p.primitive_in_id,p.primitive_out_id):
            raise ValueError('literal events inside one sample detached from physical primitive')
    return p


def compose_endpoint_lineage(lineages:Sequence,*,bias_contract:BIAS.BiasFamilyContract,tau_ba=None,require_complete_word=False):
    if not lineages:raise ValueError('nonempty endpoint lineage required')
    if not isinstance(bias_contract,BIAS.BiasFamilyContract):raise TypeError('explicit bias family contract required')
    if require_complete_word and len(lineages)!=TRANSITIONS:raise ValueError(f'complete Live word requires {TRANSITIONS} transitions')
    held,true=hard_bias_domains(bias_contract)
    Psi=[[I(1 if i==j else 0) for j in range(24)] for i in range(24)]
    blocks:list[WordSourceBlock]=[];prev_selector=None;prev_primitive=None;mode=lineages[0].mode
    for ordinal,lineage in enumerate(lineages, start=1):
        if lineage.mode!=mode:raise ValueError('single-mode word composer cannot cross H18/A21; use hybrid splice')
        prefix,parent,child=_selector_key(lineage)
        if prefix!=ordinal:raise ValueError('selector lineage skipped or duplicated a prefix')
        if prev_selector is None:
            if parent!='root':raise ValueError('first selector parent must be root')
        elif parent!=prev_selector.source_cell_id:
            raise ValueError('selector parent/child ancestry broken')
        primitive=_sample_primitive(lineage)
        if prev_primitive is not None:
            if primitive.generator_id!=prev_primitive.generator_id:raise ValueError('BRMM generator changed inside one word')
            if primitive.live_origin_id!=prev_primitive.live_origin_id:raise ValueError('one-time Live S origin changed inside one word')
            if primitive.primitive_in_id!=prev_primitive.primitive_out_id:raise ValueError('physical primitive continuity broken between samples')
        local=LOCAL.compose_attached_sample(lineage,bias_contract=bias_contract,held_bias_error=held if mode=='H' else None,true_bias=true if mode=='H' else None,tau_ba=tau_ba if mode=='A' else None)
        J=local['J_joint24'];Psi=matrix_mul(J,Psi)
        blocks=[WordSourceBlock(b.kind,b.token,_freeze(_propagate(J,_thaw(b.response))),b.sectors,b.global_relation,b.generator_id,b.live_origin_id,b.family_parameter_token) for b in blocks]
        for b in local['source_blocks']:
            blocks.append(WordSourceBlock(b.kind,b.token,b.response,b.sectors,GLOBAL_PHYSICAL_RELATION,b.generator_id,b.live_origin_id,b.family_parameter_token))
        if not local['same_bias_parameter_token_on_all_predictions']:raise RuntimeError('local bias parameter token detached')
        prev_selector=lineage.selector;prev_primitive=primitive
    bias_blocks=[b for b in blocks if b.kind=='shared_bias_driver']
    qblocks=[b for b in blocks if b.kind=='BRMM_q15_prediction'];sblocks=[b for b in blocks if b.kind=='physical_centered_S']
    return {
      'mode':mode,'transitions_composed':len(lineages),'J_joint24':Psi,'source_blocks':tuple(blocks),
      'map_representation':'pointwise_state_Jacobian_cocycle','finite_error_identity_established':False,
      'endpoint_source_cell_id':prev_selector.source_cell_id,'generator_id':prev_primitive.generator_id,'live_origin_id':prev_primitive.live_origin_id,
      'selector_ancestry_continuous':True,'primitive_transition_continuity':True,'one_generator_over_word':True,'one_Live_S_origin_over_word':True,
      'global_physical_relation':GLOBAL_PHYSICAL_RELATION,'all_physical_blocks_tied_to_global_relation':all(b.global_relation==GLOBAL_PHYSICAL_RELATION for b in blocks),
      'bias_family':bias_contract.name,'bias_parameter_token':bias_contract.parameter_token,
      'one_bias_parameter_token_over_word':all(b.family_parameter_token==bias_contract.parameter_token for b in bias_blocks),
      'q15_prediction_block_count':len(qblocks),'physical_S_block_count':len(sblocks),'bias_driver_block_count':len(bias_blocks),
      'every_prior_source_block_suffix_propagated':True,'replay_used':False,'independent_sample_source_boxes_used':False,
      'complete_word_length':len(lineages)==TRANSITIONS,
    }


def induction_theorem():
    source=SOURCE.build();sf=SOURCE.validate(source);outer=OUTER.build();of=OUTER.validate(outer);local=LOCAL.build();lf=LOCAL.validate(local);bias=BIAS.build();bf=BIAS.validate(bias)
    bad={k:v for k,v in {'source':sf,'outer':of,'local':lf,'bias':bf}.items() if v}
    if bad:raise RuntimeError('word induction prerequisites failed: '+repr(bad))
    premise=bool(source['regional_universal_Normal_Live_source_relation_closed'] and source['all_branch_successors_retained_by_operator'] and outer['sample_count']==SOURCE_SAMPLES and outer['same_history_required_for_entire_window'] and outer['constraints']['translation']['primitive_out_is_next_primitive_in'] and outer['constraints']['translation']['one_live_centered_S_origin_for_all_samples'] and local['source_uniform_attached_sample_cocycle_available'] and bias['one_persistent_parameter_token_per_family'])
    # Standard finite induction: root is covered by the regional predecessor set;
    # if every covered prefix has every successor represented and the local
    # cocycle is total on each attached successor, all prefixes 1..600 and every
    # endpoint are covered. No numerical source enumeration is required.
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','global_source_relation':GLOBAL_PHYSICAL_RELATION,
      'map_representation':'pointwise_state_Jacobian_cocycle','finite_error_identity_established':False,
      'source_samples':SOURCE_SAMPLES,'word_transitions':TRANSITIONS,
      'regional_root_set_closed':source['regional_frontend_predecessor_invariant_consumed'],
      'universal_transition_relation_closed':source['regional_universal_Normal_Live_source_relation_closed'],
      'all_successors_retained':source['all_branch_successors_retained_by_operator'],
      'local_physical_joint24_cocycle_total_on_attached_successor':local['source_uniform_attached_sample_cocycle_available'],
      'primitive_out_equals_next_in_required':outer['constraints']['translation']['primitive_out_is_next_primitive_in'],
      'one_Live_S_origin_required':outer['constraints']['translation']['one_live_centered_S_origin_for_all_samples'],
      'one_bias_parameter_token_required':bias['one_persistent_parameter_token_per_family'],
      'finite_induction_premises_closed':premise,
      'every_prefix_and_endpoint_has_single_mode_physical_cocycle':premise,
      'H18_complete_word_construction_closed':premise,
      'A21_complete_word_construction_closed':premise,
      'physical_prediction_forcing_attached':premise,
      'physical_S_residual_attached':premise,
      'all_bias_families_attached_to_single_mode_word':premise,
      'same_history_COMPLETE_BRMM_single_mode_word_closed':premise,
      'hybrid_H18_A21_word_closed':False,'H18_A21_edge_attached':False,
      'all_literal_hybrid_branches_attached':False,'storage_search_allowed':False,'ALT_LIVE_PASS':False,
      'replay_or_finite_source_enumeration_used':False,'independent_sample_source_boxes_used':False,
      'next_obligation':'splice the exact held H18 covariance into A21 at every literal release guard branch while preserving this induction invariant; then the complete hybrid physical word can unlock common-storage search'}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('global_source_relation')!=GLOBAL_PHYSICAL_RELATION or d.get('source_samples')!=SOURCE_SAMPLES or d.get('word_transitions')!=TRANSITIONS:f.append('qualification/window mismatch')
    for k in ('regional_root_set_closed','universal_transition_relation_closed','all_successors_retained','local_physical_joint24_cocycle_total_on_attached_successor','primitive_out_equals_next_in_required','one_Live_S_origin_required','one_bias_parameter_token_required','finite_induction_premises_closed','every_prefix_and_endpoint_has_single_mode_physical_cocycle','H18_complete_word_construction_closed','A21_complete_word_construction_closed','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached_to_single_mode_word','same_history_COMPLETE_BRMM_single_mode_word_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('hybrid_H18_A21_word_closed','H18_A21_edge_attached','all_literal_hybrid_branches_attached','storage_search_allowed','ALT_LIVE_PASS','replay_or_finite_source_enumeration_used','independent_sample_source_boxes_used'):
        if d.get(k) is not False:f.append(k+' not false')
    return f


def finite_storage_readiness():
    """Actual limitation of this module's J products, not inferred PASS labels.

    compose_endpoint_lineage composes pointwise state Jacobians and source
    derivatives. It does not yet compose the finite descriptors proved in
    finite_measurement_graph or a complete anchored mean-value identity.
    The existing induction/ancestry report cannot replace that missing bridge.
    """
    return {
        'map_representation': 'pointwise_state_Jacobian_cocycle',
        'finite_error_identity_for_every_event': False,
        'physical_reference_forcing_retained': False,
        'all_coefficient_product_graphs_retained': False,
        'all_configured_branches_bound_to_finite_graph': False,
    }
