#!/usr/bin/env python3
"""Outward certificate for the ALT projected common-storage inequality.

Given one interval family enclosure A in R^{n x n}, one point candidate M>0,
rho in (0,1), and an admitted coordinate supply selector C, certify

  Z' (A' M A - rho M) Z < 0   for every A in the interval family,

where Z injects ker(C).  Combined with common_storage_master's coordinate
Finsler lemma, this proves existence of a finite bounded-coordinate supply
completion for the whole enclosed family.  The arithmetic here is entirely the
repository's outward Interval + interval-LDLT backend.

This module does NOT manufacture the source-uniform endpoint family.  A caller
must provide an interval enclosure already proved to contain every physical
word under consideration.  In particular, replay matrices and finite source
samples cannot be passed off as a theorem family.
"""
from __future__ import annotations

import math
from typing import Sequence
import numpy as np

from ou3_interval import (
    Interval, matrix_point, matrix_mul, matrix_transpose,
    matrix_sub, symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
from tools.stability.ou3_alt_contraction import common_storage_master as MASTER

QUALIFICATION="OU3_ALT_PROJECTED_COMMON_STORAGE_OUTWARD_LDLT_V1"


def _shape(A):
    r=len(A); c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A): raise ValueError("ragged interval matrix")
    return r,c


def _neg(A): return [[-x for x in row] for row in A]


def _point_scaled_principal(M: np.ndarray, keep, rho: float):
    return [[Interval.outward_bounds(rho*float(M[i,j]),rho*float(M[i,j])) for j in keep] for i in keep]


def _interval_columns(A, keep):
    return [[row[j] for j in keep] for row in A]


def certify_projected_interval_family(
    A_family: Sequence[Sequence[Interval]],
    M,
    rho: float,
    supplied_indices,
) -> dict:
    """Certify the projected storage LMI for every A in one interval family."""
    n,m=_shape(A_family)
    if n==0 or n!=m: raise ValueError("A_family must be nonempty square")
    if any(not isinstance(x,Interval) for row in A_family for x in row):
        raise TypeError("A_family must be outward Interval entries")
    M=np.asarray(M,dtype=float)
    if M.shape!=(n,n) or not np.isfinite(M).all() or not np.allclose(M,M.T,rtol=0,atol=1e-12):
        raise ValueError("M must be finite symmetric n x n")
    if not math.isfinite(rho) or not 0.0<rho<1.0: raise ValueError("rho must lie in (0,1)")
    idx=tuple(int(i) for i in supplied_indices)
    if len(set(idx))!=len(idx) or any(i<0 or i>=n for i in idx): raise ValueError("invalid supplied indices")
    keep=[i for i in range(n) if i not in set(idx)]
    if not keep: raise ValueError("at least one unsupplied coordinate required")

    Mint=matrix_point(M.tolist())
    M_ok,M_piv=symmetric_positive_definite_ldlt(matrix_symmetric_hull(Mint))
    if not M_ok:
        return {"metric_spd":False,"projected_strict":False,"metric_pivots_lower":[float(x.lo) for x in M_piv],"projected_pivots_lower":[]}

    AZ=_interval_columns(A_family,keep)
    pull=matrix_mul(matrix_mul(matrix_transpose(AZ),Mint),AZ)
    rhoM=_point_scaled_principal(M,keep,rho)
    R=matrix_symmetric_hull(matrix_sub(pull,rhoM))
    ok,piv=symmetric_positive_definite_ldlt(matrix_symmetric_hull(_neg(R)))
    return {
        "metric_spd":True,
        "projected_strict":bool(ok),
        "dimension":n,
        "unsupplied_dimension":len(keep),
        "supplied_indices":list(idx),
        "rho":float(rho),
        "metric_pivots_lower":[float(x.lo) for x in M_piv],
        "projected_pivots_lower":[float(x.lo) for x in piv],
        "worst_projected_ldlt_pivot_lower":min((float(x.lo) for x in piv),default=math.inf),
        "interval_family_consumed":True,
        "ordinary_eigenvalue_used_for_certificate":False,
    }


def certify_joint24_bias_supply_family(A_family, M, rho: float) -> dict:
    return certify_projected_interval_family(A_family,M,rho,MASTER.NEUTRAL_INDICES)


def build() -> dict:
    master=MASTER.build(); mf=MASTER.validate(master)
    if mf: raise RuntimeError("neutral-supply reduction prerequisite failed: "+repr(mf))
    return {
        "qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "phase1_storage_search_allowed_consumed":master["phase1_storage_search_allowed_consumed"],
        "joint24_bias_reference_supply_selector_consumed":master["neutral_indices"],
        "outward_interval_family_arithmetic":True,
        "outward_full_matrix_LDLT_terminal_gate":True,
        "common_metric_point_candidate_only":True,
        "source_uniform_endpoint_family_must_be_proved_by_caller":True,
        "replay_or_finite_sample_family_is_qualification":False,
        "source_uniform_endpoint_family_enclosure_closed":False,
        "common_M_source_uniform_projected_LDLT_closed":False,
        "ALT_LIVE_PASS":False,
        "next_obligation":"construct a proved interval/structured enclosure of every H18 and A21 endpoint joint24 homogeneous map from the universal physical word relation, then search one M,rho and certify both families with certify_joint24_bias_supply_family",
    }


def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION or d.get("canonical_source")!="COMPLETE_BRMM_NORMAL_LIVE_WORD": f.append("qualification/source mismatch")
    for k in ("phase1_storage_search_allowed_consumed","outward_interval_family_arithmetic","outward_full_matrix_LDLT_terminal_gate","common_metric_point_candidate_only","source_uniform_endpoint_family_must_be_proved_by_caller"):
        if d.get(k) is not True: f.append(k+" not true")
    for k in ("replay_or_finite_sample_family_is_qualification","source_uniform_endpoint_family_enclosure_closed","common_M_source_uniform_projected_LDLT_closed","ALT_LIVE_PASS"):
        if d.get(k) is not False: f.append(k+" not false")
    if d.get("joint24_bias_reference_supply_selector_consumed")!=list(MASTER.NEUTRAL_INDICES): f.append("neutral selector changed")
    return f
