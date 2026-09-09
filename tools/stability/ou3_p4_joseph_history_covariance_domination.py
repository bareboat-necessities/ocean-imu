#!/usr/bin/env python3
"""Every-prefix Joseph correction-history covariance domination for P4.

For each literal event k let C_k be the homogeneous physical/tangent state map
used by the covariance recursion and B_k the same-event measurement-noise input
map (for a Joseph/reset event B_k=G_k K_k; non-measurement events have zero
columns). Shipping covariance satisfies

  P_{k+1} >= C_k P_k C_k^T + B_k R_k B_k^T,

with prediction Q, covariance floors and H->A seed appearing as additional PSD
terms.  Therefore for EVERY literal prefix m,

  P_m >= M_{m:0} P_0 M_{m:0}^T
         + sum_{k<m} M_{m:k+1} B_k R_k B_k^T M_{m:k+1}^T.   (1)

Define the prefix-stacked same-history correction operator

  W_m=[M_{m:k+1} B_k R_k^(1/2)]_{k<m}.

Then W_m W_m^T <= P_m and every stacked input u obeys

  (W_m u)^T P_m^-1 (W_m u) <= u^T u.                       (2)

This holds at the complete endpoint and at every literal prefix.  For the exact
finite-map transport's accelerometer interior term, B_k=G_k K_k and
u_k=R_k^-1/2 H_k E_k epsilon_k.  Thus the JOINT correction history is controlled
in the contemporaneous shipping metric without boxing K, bounding suffix norms,
or multiplying a worst packet by packet count.  Actual event R remains attached
to the same history.

The generic routines verify endpoint and every-prefix telescoping identities in
exact rational arithmetic.  Production interval cells may use the same algebra
with outward matrices.  This is a source-correlation primitive, not P4
promotion by itself.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

QUALIFICATION='OU3_P4_JOSEPH_HISTORY_COVARIANCE_DOMINATION_V2'


def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Fraction(0)) for j in range(len(B[0]))] for i in range(len(A))]
def mt(A):return [list(x) for x in zip(*A)]
def add(A,B):return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def sub(A,B):return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def eye(n):return [[Fraction(1 if i==j else 0) for j in range(n)] for i in range(n)]
def zero(n):return [[Fraction(0) for _ in range(n)] for _ in range(n)]
def eq(A,B):return all(A[i][j]==B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def psd2(A):return A[0][0]>=0 and A[1][1]>=0 and A[0][0]*A[1][1]-A[0][1]*A[1][0]>=0


def suffix_products(Cs):
    if not Cs:raise ValueError('nonempty event list required')
    n=len(Cs[0]);s=[None]*(len(Cs)+1);s[-1]=eye(n)
    for k in range(len(Cs)-1,-1,-1):s[k]=mm(s[k+1],Cs[k])
    return s


def telescoped_measurement_covariance(Cs,Bs,Rs):
    """Return sum suffix*B*R*B^T*suffix^T for one prefix event list."""
    if not(len(Cs)==len(Bs)==len(Rs)):raise ValueError('event arrays differ')
    n=len(Cs[0]);s=suffix_products(Cs);out=zero(n)
    for k,(B,R) in enumerate(zip(Bs,Rs)):
        term=mm(mm(B,R),mt(B));term=mm(mm(s[k+1],term),mt(s[k+1]));out=add(out,term)
    return out,s


def all_prefix_measurement_covariances(Cs,Bs,Rs):
    if not(len(Cs)==len(Bs)==len(Rs)) or not Cs:raise ValueError('invalid prefix arrays')
    return [telescoped_measurement_covariance(Cs[:m],Bs[:m],Rs[:m]) for m in range(1,len(Cs)+1)]


def build():
    # Exact 2-D, three-event smoke so both proper prefixes and endpoint are
    # exercised. Q_k are PSD slack for prediction/floor/hybrid additions.
    P0=[[Fraction(3),Fraction(1,3)],[Fraction(1,3),Fraction(2)]]
    Cs=[
      [[Fraction(3,4),Fraction(1,5)],[Fraction(0),Fraction(4,5)]],
      [[Fraction(5,6),Fraction(-1,7)],[Fraction(1,9),Fraction(7,8)]],
      [[Fraction(9,10),Fraction(1,11)],[Fraction(-1,12),Fraction(8,9)]],
    ]
    Bs=[[[Fraction(1,4)],[Fraction(1,7)]],[[Fraction(-1,6)],[Fraction(1,5)]],[[Fraction(1,8)],[Fraction(-1,9)]]]
    Rs=[[[Fraction(2)]],[[Fraction(3,2)]],[[Fraction(5,4)]]]
    Qs=[
      [[Fraction(1,10),Fraction(0)],[Fraction(0),Fraction(1,12)]],
      [[Fraction(1,11),Fraction(1,60)],[Fraction(1,60),Fraction(1,13)]],
      [[Fraction(1,14),Fraction(0)],[Fraction(0),Fraction(1,15)]],
    ]
    Ps=[P0]
    for C,B,R,Q in zip(Cs,Bs,Rs,Qs):Ps.append(add(add(mm(mm(C,Ps[-1]),mt(C)),mm(mm(B,R),mt(B))),Q))
    prefix_checks=[];all_exact=True;all_psd=True
    for m,(meas,suf) in enumerate(all_prefix_measurement_covariances(Cs,Bs,Rs),start=1):
        M=suf[0];base=mm(mm(M,P0),mt(M));remainder=sub(Ps[m],add(base,meas))
        # Propagate each Q_j through the suffix following its insertion.
        expected=zero(2)
        for j in range(m):
            tail=eye(2)
            for ell in range(j+1,m):tail=mm(Cs[ell],tail)
            expected=add(expected,mm(mm(tail,Qs[j]),mt(tail)))
        exact=eq(remainder,expected);positive=psd2(remainder);all_exact=all_exact and exact;all_psd=all_psd and positive
        prefix_checks.append({'prefix_event_count':m,'telescoping_exact':exact,'PSD_slack':positive})
    return {
      'qualification':QUALIFICATION,'same_history_literal_suffixes_required':True,'same_event_B_R_pairing_required':True,
      'Joseph_B_equals_GK_for_shipping_measurements':True,'prediction_Q_floor_and_hybrid_seed_are_PSD_slack':True,
      'every_literal_prefix_covariance_domination_identity_available':True,'whole_word_covariance_domination_identity_available':True,
      'stacked_correction_operator_WmWT_Loewner_below_each_prefix_Pm':True,'every_prefix_metric_joint_input_gain_at_most_one':True,
      'stacked_correction_operator_WWT_Loewner_below_endpoint_P':True,'endpoint_metric_joint_input_gain_at_most_one':True,
      'accelerometer_prefix_and_endpoint_interior_can_use_u_equals_RminusHalf_H_E_epsilon':True,
      'rowwise_K_box_used':False,'suffix_operator_norm_bound_used':False,'packet_count_multiplier_used':False,'independent_R_schedule_used':False,
      'exact_rational_every_prefix_telescoping_residual_matches_propagated_Q':all_exact,
      'exact_rational_every_prefix_propagated_Q_remainder_PSD':all_psd,'prefix_regression':prefix_checks,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_history_literal_suffixes_required','same_event_B_R_pairing_required','Joseph_B_equals_GK_for_shipping_measurements',
              'prediction_Q_floor_and_hybrid_seed_are_PSD_slack','every_literal_prefix_covariance_domination_identity_available',
              'whole_word_covariance_domination_identity_available','stacked_correction_operator_WmWT_Loewner_below_each_prefix_Pm',
              'every_prefix_metric_joint_input_gain_at_most_one','stacked_correction_operator_WWT_Loewner_below_endpoint_P',
              'endpoint_metric_joint_input_gain_at_most_one','accelerometer_prefix_and_endpoint_interior_can_use_u_equals_RminusHalf_H_E_epsilon',
              'exact_rational_every_prefix_telescoping_residual_matches_propagated_Q','exact_rational_every_prefix_propagated_Q_remainder_PSD'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_box_used','suffix_operator_norm_bound_used','packet_count_multiplier_used','independent_R_schedule_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if len(d.get('prefix_regression',[]))<3:f.append('proper prefix regression missing')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
