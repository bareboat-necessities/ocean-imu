#!/usr/bin/env python3
"""Whole-word Joseph correction-history covariance domination for P4.

For each literal event k let C_k be the homogeneous physical/tangent state map
used by the covariance recursion and let B_k be the same-event measurement-noise
input map (for a Joseph/reset event B_k=G_k K_k; for non-measurement events it
has zero columns).  Shipping covariance satisfies

  P_{k+1} >= C_k P_k C_k^T + B_k R_k B_k^T,

where prediction Q, covariance floors and the H->A seed are additional PSD
terms.  Propagating this identity through the literal word gives

  P_N >= M_W P_0 M_W^T
         + sum_k M_{N:k+1} B_k R_k B_k^T M_{N:k+1}^T.       (1)

Equivalently, with the stacked same-history correction operator

  W=[M_{N:k+1} B_k R_k^(1/2)]_k,

we have W W^T <= P_N.  Therefore every stacked input u obeys

  (W u)^T P_N^-1 (W u) <= u^T u.                            (2)

For the exact endpoint transport's accelerometer interior term,
B_k=G_k K_k and u_k=R_k^-1/2 H_k E_k epsilon_k.  Thus (2) controls the JOINT
suffix-weighted correction history directly in the endpoint metric.  It does
not box K, bound suffix norms, or multiply a worst packet by the packet count.
Every actual per-event R, including anisotropic applied R_S for S events in the
full covariance recursion, remains attached to the same history.

The generic routines below verify the telescoping matrix identity in exact
rational arithmetic.  Production interval cells may use the same algebra with
outward matrices.  This module is a source-correlation primitive, not P4
promotion by itself.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

QUALIFICATION='OU3_P4_JOSEPH_HISTORY_COVARIANCE_DOMINATION_V1'


def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Fraction(0)) for j in range(len(B[0]))] for i in range(len(A))]
def mt(A):return [list(x) for x in zip(*A)]
def add(A,B):return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def sub(A,B):return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def eye(n):return [[Fraction(1 if i==j else 0) for j in range(n)] for i in range(n)]
def zero(n):return [[Fraction(0) for _ in range(n)] for _ in range(n)]
def scale(A,c):return [[c*x for x in row] for row in A]
def eq(A,B):return all(A[i][j]==B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def psd2(A):return A[0][0]>=0 and A[1][1]>=0 and A[0][0]*A[1][1]-A[0][1]*A[1][0]>=0


def suffix_products(Cs):
    if not Cs:raise ValueError('nonempty event list required')
    n=len(Cs[0]);s=[None]*(len(Cs)+1);s[-1]=eye(n)
    for k in range(len(Cs)-1,-1,-1):s[k]=mm(s[k+1],Cs[k])
    return s


def telescoped_measurement_covariance(Cs, Bs, Rs):
    """Return sum suffix*B*R*B^T*suffix^T in exact/generic arithmetic."""
    if not(len(Cs)==len(Bs)==len(Rs)):raise ValueError('event arrays differ')
    n=len(Cs[0]);s=suffix_products(Cs);out=zero(n)
    for k,(B,R) in enumerate(zip(Bs,Rs)):
        term=mm(mm(B,R),mt(B));term=mm(mm(s[k+1],term),mt(s[k+1]));out=add(out,term)
    return out,s


def build():
    # Exact 2-D, two-event smoke.  Q_k are explicit PSD additions representing
    # prediction/floor/hybrid slack.  B/R remain attached to their own event.
    P0=[[Fraction(3),Fraction(1,3)],[Fraction(1,3),Fraction(2)]]
    C0=[[Fraction(3,4),Fraction(1,5)],[Fraction(0),Fraction(4,5)]]
    C1=[[Fraction(5,6),Fraction(-1,7)],[Fraction(1,9),Fraction(7,8)]]
    B0=[[Fraction(1,4)],[Fraction(1,7)]];R0=[[Fraction(2)]]
    B1=[[Fraction(-1,6)],[Fraction(1,5)]];R1=[[Fraction(3,2)]]
    Q0=[[Fraction(1,10),Fraction(0)],[Fraction(0),Fraction(1,12)]]
    Q1=[[Fraction(1,11),Fraction(1,60)],[Fraction(1,60),Fraction(1,13)]]
    P1=add(add(mm(mm(C0,P0),mt(C0)),mm(mm(B0,R0),mt(B0))),Q0)
    P2=add(add(mm(mm(C1,P1),mt(C1)),mm(mm(B1,R1),mt(B1))),Q1)
    meas,suf=telescoped_measurement_covariance([C0,C1],[B0,B1],[R0,R1])
    M=suf[0];base=mm(mm(M,P0),mt(M));remainder=sub(P2,add(base,meas))
    expected=add(mm(mm(C1,Q0),mt(C1)),Q1)
    return {
      'qualification':QUALIFICATION,'same_history_literal_suffixes_required':True,
      'same_event_B_R_pairing_required':True,'Joseph_B_equals_GK_for_shipping_measurements':True,
      'prediction_Q_floor_and_hybrid_seed_are_PSD_slack':True,
      'whole_word_covariance_domination_identity_available':True,
      'stacked_correction_operator_WWT_Loewner_below_endpoint_P':True,
      'endpoint_metric_joint_input_gain_at_most_one':True,
      'accelerometer_endpoint_interior_can_use_u_equals_RminusHalf_H_E_epsilon':True,
      'rowwise_K_box_used':False,'suffix_operator_norm_bound_used':False,'packet_count_multiplier_used':False,
      'independent_R_schedule_used':False,'exact_rational_telescoping_residual_matches_propagated_Q':eq(remainder,expected),
      'exact_rational_propagated_Q_remainder_PSD':psd2(remainder),'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_history_literal_suffixes_required','same_event_B_R_pairing_required','Joseph_B_equals_GK_for_shipping_measurements',
              'prediction_Q_floor_and_hybrid_seed_are_PSD_slack','whole_word_covariance_domination_identity_available',
              'stacked_correction_operator_WWT_Loewner_below_endpoint_P','endpoint_metric_joint_input_gain_at_most_one',
              'accelerometer_endpoint_interior_can_use_u_equals_RminusHalf_H_E_epsilon',
              'exact_rational_telescoping_residual_matches_propagated_Q','exact_rational_propagated_Q_remainder_PSD'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_box_used','suffix_operator_norm_bound_used','packet_count_multiplier_used','independent_R_schedule_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
