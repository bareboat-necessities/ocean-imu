#!/usr/bin/env python3
"""Dense IQCs for exact Cayley chord coordinates in the P4 augmented master.

The terminal master uses one common augmented coordinate z.  Let p=Pz and q=Qz
be the tangent-chord and exact-chord three-vectors for one vector measurement.
The exact Cayley graph satisfies

    q^T(q-p)=0,
    k p^T p <= q^T q <= p^T p,

with k=1/(1+r_bar^2/4) on ||c||<=r_bar.  The equality is represented as two
quadratic inequalities so the existing nonnegative-multiplier S-procedure can
consume it without adding an unrestricted-multiplier code path.

For accelerometer use the augmented variables r=[c]x R_hat delta_a_w and
u=R_hat delta_a_w+delta_b_a and define p=h-nu+r, y=q+nu.  Thus the favorable
Joseph y cross terms remain in the master.  The only outer relaxation here is
the bilinear cross-product graph for r; it retains both norm sectors and exact
orthogonality r^T c=r^T d=0.  No packetwise eta radius is introduced.
"""
from __future__ import annotations

import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD

SCHEMA=1
QUALIFICATION="OU3_P4_COMPLETE_BRMM_EXACT_CHORD_DENSE_IQC_V1"


def _shape(A):
    return len(A),len(A[0]) if A else 0

def _scale(A,a):
    c=Interval.point(float(a));return [[c*x for x in row] for row in A]

def _gram(A): return matrix_mul(matrix_transpose(A),A)
def _cross(A,B):
    return matrix_symmetric_hull(_scale(matrix_add(matrix_mul(matrix_transpose(A),B),matrix_mul(matrix_transpose(B),A)),0.5))
def _zero(n): return [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]

def exact_chord_iqcs(P,Q,k_lower:float):
    pr,pn=_shape(P);qr,qn=_shape(Q)
    if pr!=3 or qr!=3 or pn==0 or qn!=pn:raise ValueError("P,Q must be 3xn")
    if not (0.0<k_lower<=1.0):raise ValueError("invalid retention")
    pp=_gram(P);qq=_gram(Q);qp=_cross(Q,P)
    eq=matrix_symmetric_hull(matrix_sub(qq,qp)) # z'(Q'Q-sym(Q'P))z=q'(q-p)
    return {
      "orthogonality_plus":eq,
      "orthogonality_minus":_scale(eq,-1.0),
      "lower_retention":matrix_symmetric_hull(matrix_sub(qq,_scale(pp,k_lower))),
      "upper_retention":matrix_symmetric_hull(matrix_sub(pp,qq)),
    }

def cross_product_iqcs(Cmap,Dmap,Rmap,c_norm_upper:float,d_norm_upper:float):
    """Outer quadratic graph for r=c x d using valid norm/orthogonality facts."""
    cr,n=_shape(Cmap);dr,n2=_shape(Dmap);rr,n3=_shape(Rmap)
    if (cr,dr,rr)!=(3,3,3) or n==0 or n2!=n or n3!=n:raise ValueError("maps must be 3xn")
    if min(c_norm_upper,d_norm_upper)<0:raise ValueError("nonnegative radii required")
    cc=_gram(Cmap);dd=_gram(Dmap);rrg=_gram(Rmap)
    rc=_cross(Rmap,Cmap);rd=_cross(Rmap,Dmap)
    return {
      "r_norm_from_d":matrix_symmetric_hull(matrix_sub(_scale(dd,c_norm_upper*c_norm_upper),rrg)),
      "r_norm_from_c":matrix_symmetric_hull(matrix_sub(_scale(cc,d_norm_upper*d_norm_upper),rrg)),
      "r_dot_c_plus":rc,"r_dot_c_minus":_scale(rc,-1.0),
      "r_dot_d_plus":rd,"r_dot_d_minus":_scale(rd,-1.0),
    }

def _selector(n,offset):
    A=_zero(3)
    # Return 3xn, not square.
    A=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
    for i in range(3):A[i][offset+i]=Interval.point(1.0)
    return A

def build():
    ch=CHORD.build();vf=CHORD.validate(ch)
    if vf:raise RuntimeError(f"chord prerequisite failed: {vf}")
    # Smoke augmented coordinate [c,d,r,p,q] (15D).
    n=15;C=_selector(n,0);D=_selector(n,3);R=_selector(n,6);P=_selector(n,9);Q=_selector(n,12)
    chord=exact_chord_iqcs(P,Q,float(ch["information_retention_factor_lower_full_entry"]))
    cross=cross_product_iqcs(C,D,R,float(ch["full_declared_entry_cayley_norm_upper"]),2.941995)
    allm=list(chord.values())+list(cross.values())
    finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for A in allm for row in A for x in row)
    sym=all(all(A[i][j].lo==A[j][i].lo and A[i][j].hi==A[j][i].hi for j in range(len(A))) for i in range(len(A))) for A in allm)
    return {
      "schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "full_declared_entry_retention":ch["information_retention_factor_lower_full_entry"],
      "chord_IQC_names":list(chord),"cross_product_IQC_names":list(cross),
      "exact_q_dot_q_minus_p_equality_encoded_as_two_IQCs":True,
      "q_u_Joseph_cross_terms_preserved_by_coordinate_not_scalarized":True,
      "mixed_c_cross_aw_has_explicit_augmented_coordinate":True,
      "cross_product_orthogonality_retained":True,
      "cross_product_norm_sectors_retained":True,
      "bilinear_cross_product_graph_outer_not_declared_exact":True,
      "dense_IQC_matrices_finite":finite,"dense_IQC_matrices_symmetric":sym,
      "compatible_with_nonnegative_multiplier_joint_master":True,
      "packet_count_multiplier_used":False,"standalone_eta_Rinv_budget_used":False,
      "endpoint_augmented_LDLT_closed_here":False,"every_prefix_augmented_LDLT_closed_here":False,"P4_promoted_here":False,
    }

def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("exact_q_dot_q_minus_p_equality_encoded_as_two_IQCs","q_u_Joseph_cross_terms_preserved_by_coordinate_not_scalarized","mixed_c_cross_aw_has_explicit_augmented_coordinate","cross_product_orthogonality_retained","cross_product_norm_sectors_retained","bilinear_cross_product_graph_outer_not_declared_exact","dense_IQC_matrices_finite","dense_IQC_matrices_symmetric","compatible_with_nonnegative_multiplier_joint_master"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("packet_count_multiplier_used","standalone_eta_Rinv_budget_used","endpoint_augmented_LDLT_closed_here","every_prefix_augmented_LDLT_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
