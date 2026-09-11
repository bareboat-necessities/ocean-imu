#!/usr/bin/env python3
"""Outer induction from same-cell local Jacobian families to complete-word maps.

This is the missing set-theoretic bridge between Phase-1's universal physical
word relation and the projected common-storage LDLT backend.

At literal event k let F_k be a finite family of outward interval matrices such
that every *actual same-history* local homogeneous Jacobian belongs to at least
one member of F_k.  If E_k encloses every reachable prefix product through k,
then

    E_{k+1} = hull { A E : A in F_{k+1}, E in E_k }

contains every reachable product through k+1.  The hull is an OUTER inclusion
only.  It is never interpreted as a realizable source history, and it is formed
only after each local map has been generated from one correlated estimator-owned
source cell.  Hence coefficient/source correlations are not broken in the
local shipping operation; only the image set is conservatively enlarged.

Keeping several partition labels avoids forcing unrelated modes/event geometries
into one Cartesian matrix box.  Labels have no theorem semantics beyond set
partitioning: every successor must name a parent partition and every local map
must already have a source-uniform same-cell qualification.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul
from tools.stability.ou3_alt_contraction import phase1_closure as PHASE1

QUALIFICATION="OU3_ALT_ENDPOINT_FAMILY_OUTER_INDUCTION_V1"


def _shape(A):
    r=len(A); c=len(A[0]) if r else 0
    if any(len(x)!=c for x in A): raise ValueError("ragged matrix")
    return r,c


def _check_interval_matrix(A,n=None):
    r,c=_shape(A)
    if r==0 or r!=c or (n is not None and r!=n): raise ValueError("matrix must be nonempty square")
    if any(not isinstance(v,Interval) for row in A for v in row): raise TypeError("matrix entries must be outward Interval")
    return r


def hull_interval(a:Interval,b:Interval)->Interval:
    return Interval(min(a.lo,b.lo),max(a.hi,b.hi))


def hull_matrices(mats:Sequence[Sequence[Sequence[Interval]]]):
    if not mats: raise ValueError("cannot hull empty matrix family")
    n=_check_interval_matrix(mats[0])
    for A in mats[1:]: _check_interval_matrix(A,n)
    out=[[mats[0][i][j] for j in range(n)] for i in range(n)]
    for A in mats[1:]:
        out=[[hull_interval(out[i][j],A[i][j]) for j in range(n)] for i in range(n)]
    return out


def matrix_encloses(outer,inner)->bool:
    if _shape(outer)!=_shape(inner): return False
    return all(outer[i][j].lo<=inner[i][j].lo and outer[i][j].hi>=inner[i][j].hi
               for i in range(len(outer)) for j in range(len(outer[0])))


@dataclass(frozen=True)
class QualifiedLocalFamily:
    label:str
    matrices:tuple
    source_uniform_same_cell_qualified:bool
    replay_or_finite_sample_derived:bool=False

    def validate(self,n:int)->None:
        if not isinstance(self.label,str) or not self.label: raise ValueError("local family label required")
        if not self.matrices: raise ValueError("local family cannot be empty")
        for A in self.matrices:_check_interval_matrix(A,n)
        if self.source_uniform_same_cell_qualified is not True:
            raise ValueError("local family lacks source-uniform same-cell qualification")
        if self.replay_or_finite_sample_derived:
            raise ValueError("replay/finite-sample family cannot qualify theorem induction")


@dataclass(frozen=True)
class PartitionSuccessor:
    child_label:str
    parent_labels:tuple[str,...]
    local_family:QualifiedLocalFamily


def initial_partition(n:int,label="root"):
    return {label:matrix_identity(n)}


def advance_partitions(prefixes:dict[str,object], successors:Sequence[PartitionSuccessor]):
    """One exact set-inclusion induction step.

    Each child hull contains the union of A*E for every named parent E and every
    locally-qualified A.  Missing parents and unqualified maps fail closed.
    """
    if not prefixes: raise ValueError("empty predecessor partition")
    n=None
    for E in prefixes.values(): n=_check_interval_matrix(E,n)
    if not successors: raise ValueError("empty successor relation")
    grouped:dict[str,list]= {}
    for s in successors:
        if not isinstance(s.child_label,str) or not s.child_label: raise ValueError("child label required")
        if not s.parent_labels: raise ValueError("successor must retain at least one parent")
        s.local_family.validate(n)
        images=[]
        for p in s.parent_labels:
            if p not in prefixes: raise ValueError("successor names unreachable/missing parent "+repr(p))
            E=prefixes[p]
            for A in s.local_family.matrices: images.append(matrix_mul(A,E))
        grouped.setdefault(s.child_label,[]).extend(images)
    return {label:hull_matrices(images) for label,images in grouped.items()}


def prove_endpoint_outer_induction(n:int, steps:Sequence[Sequence[PartitionSuccessor]], root_label="root"):
    """Construct endpoint outer partitions and an explicit induction ledger."""
    if n<=0: raise ValueError("positive dimension required")
    prefixes=initial_partition(n,root_label); ledger=[]
    for k,succ in enumerate(steps):
        before=tuple(sorted(prefixes))
        prefixes=advance_partitions(prefixes,succ)
        after=tuple(sorted(prefixes))
        ledger.append({"step":k,"parents":before,"children":after,
                       "all_successors_retained":True,"outer_hull_is_realizable_history":False})
    return {"endpoint_partitions":prefixes,"ledger":ledger,
            "steps":len(steps),"outer_union_inclusion_by_induction":True,
            "outer_hull_is_realizable_source_history":False}


def build():
    p=PHASE1.build(); f=PHASE1.validate(p)
    if f: raise RuntimeError("Phase-1 prerequisite failed: "+repr(f))
    return {
      "qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "phase1_storage_search_allowed_consumed":bool(p["storage_search_allowed"]),
      "local_maps_must_be_estimator_owned_same_cell_before_hull":True,
      "replay_or_finite_sample_local_family_forbidden":True,
      "entrywise_hull_used_only_as_outer_union_inclusion":True,
      "outer_hull_may_generate_theorem_history":False,
      "all_successors_required_each_induction_step":True,
      "partitioned_hulls_supported":True,
      "endpoint_outer_inclusion_induction_materialized":True,
      "actual_600_step_source_family_bound":False,
      "ALT_LIVE_PASS":False,
      "next_obligation":"instantiate every step's QualifiedLocalFamily from the universal estimator-owned source-cell relation for H18 and A21, retaining all literal guard/event partitions; then feed each endpoint partition to projected_storage_certificate with one common M,rho",
    }


def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION or d.get("canonical_source")!="COMPLETE_BRMM_NORMAL_LIVE_WORD":f.append("qualification/source mismatch")
    for k in ("phase1_storage_search_allowed_consumed","local_maps_must_be_estimator_owned_same_cell_before_hull","replay_or_finite_sample_local_family_forbidden","entrywise_hull_used_only_as_outer_union_inclusion","all_successors_required_each_induction_step","partitioned_hulls_supported","endpoint_outer_inclusion_induction_materialized"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("outer_hull_may_generate_theorem_history","actual_600_step_source_family_bound","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return f
