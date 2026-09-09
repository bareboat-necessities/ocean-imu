#!/usr/bin/env python3
"""Exact finite-angle chord coordinate for the complete-BRMM signed Joseph ledger.

This module deliberately keeps the *vector* chord coordinate, not only a norm
bound.  For Cayley c, C=[c]x, den=1+||c||^2/4 and E the exact rotation,

    (E-I)m = (I + C/2) (C m) / den.                         (1)

Hence

    ||(E-I)m||^2 = ||C m||^2 / den.                         (2)

For the shipping accelerometer let d=R_hat*delta_a_w and b=delta_b_a.  The
exact physical residual is

    y_a = (E-I) f_hat + E d + b
        = (E-I)(f_hat+d) + d + b.                           (3)

Define p=C(f_hat+d), q=(E-I)(f_hat+d), u=d+b.  Then y_a=q+u and
q=(I+C/2)p/den exactly.  The bilinear term C d is therefore retained inside p;
it is NOT charged as an independent eta radius.  Because q itself is retained,
the Joseph information keeps the q-u cross term in y_a^T S^-1 y_a.

The same construction applies to magnetometer with u=0.  The information
retention relative to ||C m||^2 is 1/den = cos(theta/2)^2.  It must not be
squared again.

This is an algebra/graph-coordinate producer.  Source-uniform P/H/R/K cells,
BIAS1 continuation and reset transport are attached by downstream same-history
lineage consumers.  It never promotes P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_add, matrix_identity, matrix_mul
import ou3_p4_cayley_sector_certificate as CAYLEY

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_EXACT_CHORD_JOINT_COORDINATE_V1"


def _skew(c: Sequence[Sequence[Interval]]):
    if len(c) != 3 or any(len(row) != 1 for row in c):
        raise ValueError("c must be 3x1")
    z=Interval.point(0.0); x,y,w=(row[0] for row in c)
    return [[z,-w,y],[w,z,-x],[-y,x,z]]


def exact_rotation(c):
    C=_skew(c); C2=matrix_mul(C,C)
    r2=sum((row[0].square() for row in c), Interval.point(0.0))
    den=Interval.point(1.0)+Interval.point(0.25)*r2
    EmI=[[ (C[i][j]+Interval.point(0.5)*C2[i][j])/den for j in range(3)] for i in range(3)]
    return matrix_add(matrix_identity(3),EmI),EmI,C,den


def chord_coordinate(c,m):
    if len(m)!=3 or any(len(row)!=1 for row in m): raise ValueError("m must be 3x1")
    E,EmI,C,den=exact_rotation(c)
    p=matrix_mul(C,m)
    q=matrix_mul(EmI,m)
    # Exact vector expression q=(I+C/2)p/den.
    A=matrix_add(matrix_identity(3), [[Interval.point(0.5)*C[i][j] for j in range(3)] for i in range(3)])
    q_from_p=[[x/den for x in row] for row in matrix_mul(A,p)]
    return {"E":E,"C":C,"den":den,"p":p,"q":q,"q_from_p":q_from_p}


def accelerometer_joint_coordinate(c,f_hat,R_hat,delta_aw,delta_ba):
    for name,v in (("f_hat",f_hat),("delta_aw",delta_aw),("delta_ba",delta_ba)):
        if len(v)!=3 or any(len(row)!=1 for row in v): raise ValueError(name+" must be 3x1")
    if len(R_hat)!=3 or any(len(row)!=3 for row in R_hat): raise ValueError("R_hat must be 3x3")
    d=matrix_mul(R_hat,delta_aw)
    m=matrix_add(f_hat,d)
    ch=chord_coordinate(c,m)
    u=matrix_add(d,delta_ba)
    y=matrix_add(ch["q"],u)
    # Direct shipping identity: (E-I)f + E d + b.
    direct=matrix_add(matrix_add(matrix_mul(matrix_add(ch["E"], [[-x for x in row] for row in matrix_identity(3)]),f_hat),matrix_mul(ch["E"],d)),delta_ba)
    return {**ch,"d_body":d,"m_joint":m,"u_linear":u,"y_joint":y,"y_direct":direct}


def _contains_zero(A):
    return all(x.lo<=0.0<=x.hi for row in A for x in row)


def _sub(A,B):
    return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def build(domain_path: Path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve(); domain=json.loads(path.read_text())
    cay=CAYLEY.build(path); vf=CAYLEY.validate(cay)
    if vf: raise RuntimeError(f"Cayley prerequisite failed: {vf}")
    theta_deg=float(domain["initial_filter_entrance"]["attitude"]["full_attitude_error_upper_deg"])
    theta=math.radians(theta_deg); r=2.0*math.tan(theta/2.0)
    retention=math.nextafter(1.0/(1.0+0.25*r*r),-math.inf)
    # Non-axis point smoke exercises the exact vector and accelerometer identities.
    pt=lambda x:[[Interval.point(float(v))] for v in x]
    c=pt([0.31,-0.17,0.22]); f=pt([1.4,-0.8,-9.2]); da=pt([0.4,-0.3,0.2]); ba=pt([0.05,-0.02,0.01])
    R=[[Interval.point(1.0 if i==j else 0.0) for j in range(3)] for i in range(3)]
    a=accelerometer_joint_coordinate(c,f,R,da,ba)
    vec_ok=_contains_zero(_sub(a["q"],a["q_from_p"]))
    acc_ok=_contains_zero(_sub(a["y_joint"],a["y_direct"]))
    return {
      "schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "filter_changed":False,"declared_domain_changed":False,"trajectory_replay_used":False,
      "full_declared_entry_attitude_deg":theta_deg,"full_declared_entry_cayley_norm_upper":r,
      "exact_chord_vector_identity_closed":vec_ok,"exact_accelerometer_joint_identity_closed":acc_ok,
      "accelerometer_joint_vector":"m=f_hat+R_hat*delta_a_w",
      "accelerometer_linear_companion":"u=R_hat*delta_a_w+delta_b_a",
      "exact_residual":"y=q+u; q=(I+[c]x/2)[c]x*m/(1+||c||^2/4)",
      "mixed_c_cross_aw_retained_inside_chord_coordinate":True,
      "accelerometer_bias_retained_linearly_in_same_residual":True,
      "joseph_q_u_cross_term_must_be_retained":True,
      "standalone_eta_Rinv_budget_used":False,"packet_count_multiplier_used":False,
      "information_retention_factor_lower_full_entry":retention,
      "information_retention_is_cos_half_squared_not_fourth_power":True,
      "full_45deg_information_retention_gt_0_85":retention>0.85,
      "same_history_P_H_R_K_attachment_required_downstream":True,
      "physical_BIAS1_attachment_required_downstream":True,
      "finite_reset_transport_required_downstream":True,
      "endpoint_augmented_LDLT_closed_here":False,"every_prefix_augmented_LDLT_closed_here":False,"P4_promoted_here":False,
    }


def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("exact_chord_vector_identity_closed","exact_accelerometer_joint_identity_closed","mixed_c_cross_aw_retained_inside_chord_coordinate","accelerometer_bias_retained_linearly_in_same_residual","joseph_q_u_cross_term_must_be_retained","information_retention_is_cos_half_squared_not_fourth_power","full_45deg_information_retention_gt_0_85","same_history_P_H_R_K_attachment_required_downstream","physical_BIAS1_attachment_required_downstream","finite_reset_transport_required_downstream"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("filter_changed","declared_domain_changed","trajectory_replay_used","standalone_eta_Rinv_budget_used","packet_count_multiplier_used","endpoint_augmented_LDLT_closed_here","every_prefix_augmented_LDLT_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    if float(d.get("full_declared_entry_attitude_deg",0))!=45.0:f.append("full entry is not 45 deg")
    r=float(d.get("information_retention_factor_lower_full_entry",0))
    if not (0.85<r<1.0):f.append("invalid full-entry chord retention")
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args()
    d=build(a.domain);f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"retention":d["information_retention_factor_lower_full_entry"],"vector_identity":d["exact_chord_vector_identity_closed"],"accelerometer_identity":d["exact_accelerometer_joint_identity_closed"],"failures":f},indent=2));return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
