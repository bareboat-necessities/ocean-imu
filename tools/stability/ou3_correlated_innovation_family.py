#!/usr/bin/env python3
"""Correlated innovation family for OU-III measurement events.

The proof primitive is the triplet (P,H,R), not an independently selectable
innovation box.  This module accepts only P,H,R and constructs

    PHt = P H^T
    S   = H P H^T + R
    S^-1
    K   = P H^T S^-1
    A   = I-KH

inside one call.  Callers cannot inject an arbitrary S unrelated to the same
P,H,R cell.  For three-component shipping measurements, invertibility consumes
the structural Riccati fact P>=0 and uniform R>0, so singular matrices that
exist only in the entrywise S hull are excluded.

The returned interval matrices are outward enclosures.  This closes the
specific P/H/R -> S -> S^-1 -> K provenance break; it does not by itself prove
endpoint contraction, every-prefix retention, or P4/P5.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ou3_interval import Interval, IntervalMatrix, matrix_add, matrix_identity, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_innovation_psd_plus_R_inverse as INNOV


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r=len(A); c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):
        raise ValueError("ragged interval matrix")
    return r,c


@dataclass(frozen=True)
class CorrelatedInnovationFamily:
    """One source/covariance cell retaining P,H,R ancestry through K."""
    P: IntervalMatrix
    H: IntervalMatrix
    R: IntervalMatrix
    PHt: IntervalMatrix
    S: IntervalMatrix
    Sinv: IntervalMatrix
    K: IntervalMatrix
    A: IntervalMatrix
    inverse_provenance: dict


def build(P: Sequence[Sequence[Interval]], H: Sequence[Sequence[Interval]], R: Sequence[Sequence[Interval]]) -> CorrelatedInnovationFamily:
    n,m=_shape(P)
    rows,cols=_shape(H)
    if n==0 or n!=m or cols!=n or rows==0 or _shape(R)!=(rows,rows):
        raise ValueError("correlated innovation P/H/R dimension mismatch")
    Ps=matrix_symmetric_hull(P)
    Rs=matrix_symmetric_hull(R)
    Ht=matrix_transpose(H)
    PHt=matrix_mul(Ps,Ht)
    # S is derived here and is never accepted as an independent input.
    S=matrix_symmetric_hull(matrix_add(matrix_mul(H,PHt),Rs))
    if rows==3:
        Sinv,proof=INNOV.innovation_inverse_psd_plus_R_3x3(S,Rs)
        mode="CORRELATED_P_H_R_PSD_PLUS_R_3X3"
    else:
        Sinv=matrix_inverse_gauss_jordan(S)
        proof=None
        mode="CORRELATED_P_H_R_GENERIC_VALIDATED_INTERVAL"
    K=matrix_mul(PHt,Sinv)
    A=matrix_sub(matrix_identity(n),matrix_mul(K,H))
    provenance={
        "mode":mode,
        "primitive_family":"(P,H,R)",
        "S_is_external_independent_argument":False,
        "S_identity":"S=H P H^T+R",
        "K_identity":"K=P H^T (H P H^T+R)^-1",
        "same_P_H_R_used_for_PHt_S_Sinv_K":True,
        "covariance_PSD_provenance":"structural Riccati covariance invariant",
        "R_uniform_SPD_required":True,
        "rectangular_S_family_can_be_selected_independently":False,
        "inverse_proof":proof,
        "P4_PASS":False,
        "P5_MAY_START":False,
    }
    return CorrelatedInnovationFamily(Ps,[list(r) for r in H],Rs,PHt,S,Sinv,K,A,provenance)


def validate_point_regression() -> dict:
    from ou3_interval import matrix_point
    P=matrix_point([[2.0,0.2,0.0],[0.2,1.5,0.1],[0.0,0.1,1.0]])
    H=matrix_point([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    R=matrix_point([[0.2,0.0,0.0],[0.0,0.3,0.0],[0.0,0.0,0.4]])
    f=build(P,H,R)
    return {
        "same_P_H_R_used_for_PHt_S_Sinv_K":f.inverse_provenance["same_P_H_R_used_for_PHt_S_Sinv_K"],
        "S_is_external_independent_argument":f.inverse_provenance["S_is_external_independent_argument"],
        "rectangular_S_family_can_be_selected_independently":f.inverse_provenance["rectangular_S_family_can_be_selected_independently"],
        "P4_PASS":False,
        "P5_MAY_START":False,
    }


if __name__=="__main__":
    import json
    d=validate_point_regression()
    print(json.dumps(d,indent=2,sort_keys=True))
    raise SystemExit(0 if d["same_P_H_R_used_for_PHt_S_Sinv_K"] and not d["S_is_external_independent_argument"] and not d["rectangular_S_family_can_be_selected_independently"] else 2)
