#!/usr/bin/env python3
"""Exact same-Joseph-cell identities used by the finite-reset P4 ledger.

For SPD P,R, H, S=HPH'+R, K=PH'S^-1 and P+=P-KSK',

    P+^-1 K = H' R^-1,
    K' P+^-1 K = R^-1 - S^-1.                         (1)

For a physical residual y=H e+eta, d=Ky and t=e-d,

    t' P+^-1 d = y' S^-1 y - eta' R^-1 y.             (2)

These identities are important because they retain the actual K/y correlation
inside the signed Joseph/reset ledger.  They are matrix equalities, not bounds,
and do not require a covariance condition-number scalarization.

This module supplies exact-fraction regression checks and documents the form to
be consumed by interval event cells.  It does not promote P4 by itself.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

QUALIFICATION='OU3_P4_SAME_CELL_KALMAN_CROSS_IDENTITIES_V1'

def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Fraction(0)) for j in range(len(B[0]))] for i in range(len(A))]
def mt(A):return [list(x) for x in zip(*A)]
def mv(A,x):return [sum((A[i][j]*x[j] for j in range(len(x))),Fraction(0)) for i in range(len(A))]
def add(A,B):return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def sub(A,B):return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Fraction(0))
def inv2(A):
    d=A[0][0]*A[1][1]-A[0][1]*A[1][0]
    if d==0:raise ValueError('singular')
    return [[A[1][1]/d,-A[0][1]/d],[-A[1][0]/d,A[0][0]/d]]
def eye(n):return [[Fraction(1 if i==j else 0) for j in range(n)] for i in range(n)]
def zero(A):return all(x==0 for row in A for x in row)

def build():
    # Nonsymmetric H and dense SPD P/R smoke, exact rational arithmetic.
    P=[[Fraction(3),Fraction(1,2)],[Fraction(1,2),Fraction(2)]]
    H=[[Fraction(1),Fraction(2)],[Fraction(-1),Fraction(1)]]
    R=[[Fraction(2),Fraction(1,4)],[Fraction(1,4),Fraction(3)]]
    Pinv=inv2(P);Rinv=inv2(R);A=mm(mm(H,P),mt(H));S=add(A,R);Sinv=inv2(S);K=mm(mm(P,mt(H)),Sinv)
    Pplus=sub(P,mm(mm(K,S),mt(K)));Jplus=inv2(Pplus)
    id1=sub(mm(Jplus,K),mm(mt(H),Rinv))
    id2=sub(mm(mm(mt(K),Jplus),K),sub(Rinv,Sinv))
    e=[Fraction(2,3),Fraction(-1,5)];eta=[Fraction(1,7),Fraction(-1,11)];y=[a+b for a,b in zip(mv(H,e),eta)];d=mv(K,y);t=[a-b for a,b in zip(e,d)]
    lhs=dot(t,mv(Jplus,d));rhs=dot(y,mv(Sinv,y))-dot(eta,mv(Rinv,y));id3=lhs-rhs
    info_form=sub(Jplus,add(Pinv,mm(mm(mt(H),Rinv),H)))
    return {'qualification':QUALIFICATION,'same_cell_required':True,'Pplus_inverse_K_equals_Ht_Rinverse':zero(id1),'Kt_Pplus_inverse_K_equals_Rinverse_minus_Sinverse':zero(id2),'posterior_correction_cross_identity_closed':id3==0,'posterior_information_form_identity_closed':zero(info_form),'correction_energy_tied_to_same_Joseph_information':True,'condition_number_scalarization_used':False,'independent_K_box_used':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_required','Pplus_inverse_K_equals_Ht_Rinverse','Kt_Pplus_inverse_K_equals_Rinverse_minus_Sinverse','posterior_correction_cross_identity_closed','posterior_information_form_identity_closed','correction_energy_tied_to_same_Joseph_information'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('condition_number_scalarization_used','independent_K_box_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
