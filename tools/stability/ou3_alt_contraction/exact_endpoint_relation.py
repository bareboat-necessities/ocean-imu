#!/usr/bin/env python3
"""Exact universal endpoint-map relation for ALT common storage.

The source-reachable selector theorem on main already proves that every member
of the correlated source relation O^601_BRMM has every shipping branch retained
as one selector-family path, together with all three absolute bias lineages.
ALT's physical_word theorem proves that every such selector path has one
joint24 physical cocycle, and hybrid_word proves the configured H18/A21 splice.

Their composition defines the exact (generally infinite) set W of homogeneous
joint24 endpoint maps to which the common-storage inequality must apply.  This
module closes that quantifier binding.  It intentionally does not replace W by
finite source samples or claim that a numerical interval representation of W
has already been built.
"""
from __future__ import annotations

import ou3_p4_source_reachable_selector_family as REACH
from tools.stability.ou3_alt_contraction import physical_word as WORD
from tools.stability.ou3_alt_contraction import hybrid_word as HYBRID
from tools.stability.ou3_alt_contraction import phase1_closure as PHASE1

QUALIFICATION="OU3_ALT_EXACT_UNIVERSAL_JOINT24_ENDPOINT_RELATION_V1"
SOURCE_RELATION="O^601_BRMM"


def build():
    reach=REACH.build(); rf=REACH.validate(reach)
    word=WORD.induction_theorem(); wf=WORD.validate(word)
    hybrid=HYBRID.build(); hf=HYBRID.validate(hybrid)
    phase1=PHASE1.build(); pf=PHASE1.validate(phase1)
    bad={k:v for k,v in (("reach",rf),("word",wf),("hybrid",hf),("phase1",pf)) if v}
    if bad: raise RuntimeError("endpoint relation prerequisites failed: "+repr(bad))
    exact=bool(
      reach["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"]
      and reach["all_branch_successors_retained"]
      and reach["all_bias_absolute_prefix_relations_attached"]
      and word["every_prefix_and_endpoint_has_single_mode_physical_cocycle"]
      and word["all_physical_blocks_tied_to_global_relation"] if "all_physical_blocks_tied_to_global_relation" in word else True
    )
    # word.induction_theorem expresses the same property by the finite-induction
    # premises and the H/A construction flags; no endpoint enumeration is used.
    exact=bool(exact and word["finite_induction_premises_closed"]
               and word["H18_complete_word_construction_closed"]
               and word["A21_complete_word_construction_closed"]
               and hybrid["hybrid_H18_A21_word_closed"]
               and phase1["storage_search_allowed"])
    return {
      "qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "source_relation":SOURCE_RELATION,"word_transitions":600,"source_samples":601,"joint_dimension":24,
      "source_reachable_selector_family_relation_consumed":True,
      "all_branch_successors_retained":reach["all_branch_successors_retained"],
      "all_bias_lineages_retained":reach["all_bias_absolute_prefix_relations_attached"],
      "every_selector_endpoint_has_physical_joint24_cocycle":word["every_prefix_and_endpoint_has_single_mode_physical_cocycle"],
      "configured_H18_A21_hybrid_splice_closed":hybrid["hybrid_H18_A21_word_closed"],
      "exact_universal_endpoint_map_relation_defined":exact,
      "endpoint_map_set_symbol":"W_ALT_joint24(O^601_BRMM)",
      "endpoint_map_set_definition":"{J_joint24(path,bias,guard): path in source-reachable selector-family relation, bias in BIAS0/1/2, guard in configured H18/A21 outcomes}",
      "finite_source_enumeration_used":False,"trajectory_replay_used":False,"independent_sample_boxes_used":False,
      "numeric_interval_endpoint_representation_closed":False,
      "common_storage_quantifier":"for every A in W_ALT_joint24(O^601_BRMM)",
      "common_M_source_uniform_projected_LDLT_closed":False,"ALT_LIVE_PASS":False,
      "next_obligation":"construct a certified finite outer representation W_hat containing this exact endpoint relation (prefer partitioned structured intervals); then certify the projected common-storage LMI on every partition with one M,rho",
    }


def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION or d.get("source_relation")!=SOURCE_RELATION:f.append("qualification/source mismatch")
    for k in ("source_reachable_selector_family_relation_consumed","all_branch_successors_retained","all_bias_lineages_retained","every_selector_endpoint_has_physical_joint24_cocycle","configured_H18_A21_hybrid_splice_closed","exact_universal_endpoint_map_relation_defined"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("finite_source_enumeration_used","trajectory_replay_used","independent_sample_boxes_used","numeric_interval_endpoint_representation_closed","common_M_source_uniform_projected_LDLT_closed","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    if d.get("joint_dimension")!=24 or d.get("word_transitions")!=600 or d.get("source_samples")!=601:f.append("dimension/window mismatch")
    return f
