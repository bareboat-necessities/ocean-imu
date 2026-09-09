#!/usr/bin/env python3
"""Exact reduced same-cell Joseph+reset signed identity for P4.

Let one shipping Joseph cell have prior P, precision J=P^-1, tangent H,
measurement covariance R, S=H P H^T+R, K=P H^T S^-1, and Joseph posterior
P+.  For the exact physical residual

    q = H e + eta,
    d = K q,
    t = e-d,

and the reset defect expressed in the Joseph frame

    b = G^-1 rho,

the ordinary signed ledger is

  Delta V = -q^T S^-1 q + eta^T R^-1 eta
            +2 t^T (P+)^-1 b + b^T (P+)^-1 b.

Using the exact same-cell identities

  (P+)^-1 K = H^T R^-1,
  (P+)^-1 = J + H^T R^-1 H,

this is IDENTICALLY

  Delta V = -q^T R^-1 q + eta^T R^-1 eta
            +q^T R^-1 H d
            +2 e^T J b - 2 eta^T R^-1 H b
            +b^T J b + (H b)^T R^-1 (H b).            (*)

Thus the production augmented master does not need an independently enclosed
S^-1 or posterior precision for nonlinear/reset energy.  K survives only in the
same-cell graph coordinate d=Kq.  The prior moving precision J is the same
shipping metric already carried by the complete word.

This module verifies (*) in exact rational arithmetic.  It is an algebraic
reduction, not a P4 promotion.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

import ou3_p4_same_cell_kalman_cross_identities as X

QUALIFICATION='OU3_P4_REDUCED_SAME_CELL_JOSEPH_RESET_SIGNED_IDENTITY_V1'


def quad(A,x):return X.dot(x,X.mv(A,x))

def build():
    P=[[Fraction(3),Fraction(1,2)],[Fraction(1,2),Fraction(2)]]
    H=[[Fraction(1),Fraction(2)],[Fraction(-1),Fraction(1)]]
    R=[[Fraction(2),Fraction(1,4)],[Fraction(1,4),Fraction(3)]]
    J=X.inv2(P);Rinv=X.inv2(R);S=X.add(X.mm(X.mm(H,P),X.mt(H)),R);Sinv=X.inv2(S)
    K=X.mm(X.mm(P,X.mt(H)),Sinv);Pplus=X.sub(P,X.mm(X.mm(K,S),X.mt(K)));Jplus=X.inv2(Pplus)
    e=[Fraction(2,3),Fraction(-1,5)];eta=[Fraction(1,7),Fraction(-1,11)]
    q=[a+b for a,b in zip(X.mv(H,e),eta)];d=X.mv(K,q);t=[a-b for a,b in zip(e,d)]
    # b is already the reset defect after G^-1; choose a dense exact smoke.
    b=[Fraction(1,13),Fraction(-2,17)]
    original=-quad(Sinv,q)+quad(Rinv,eta)+2*X.dot(t,X.mv(Jplus,b))+quad(Jplus,b)
    Hb=X.mv(H,b)
    reduced=(-quad(Rinv,q)+quad(Rinv,eta)+X.dot(q,X.mv(Rinv,X.mv(H,d)))
             +2*X.dot(e,X.mv(J,b))-2*X.dot(eta,X.mv(Rinv,Hb))
             +quad(J,b)+quad(Rinv,Hb))
    joseph_reduction=(-quad(Sinv,q)+quad(Rinv,eta))-(-quad(Rinv,q)+quad(Rinv,eta)+X.dot(q,X.mv(Rinv,X.mv(H,d))))
    reset_reduction=(2*X.dot(t,X.mv(Jplus,b))+quad(Jplus,b))-(2*X.dot(e,X.mv(J,b))-2*X.dot(eta,X.mv(Rinv,Hb))+quad(J,b)+quad(Rinv,Hb))
    return {
      'qualification':QUALIFICATION,'same_shipping_P_H_R_K_cell_required':True,
      'same_cell_cross_identity_consumed':True,'S_inverse_eliminated_from_reduced_nonlinear_reset_ledger':True,
      'posterior_precision_eliminated_from_reduced_nonlinear_reset_ledger':True,
      'prior_precision_is_same_moving_shipping_metric':True,'K_survives_only_through_same_cell_d_equals_Kq_graph':True,
      'reduced_identity':'DeltaV=-q^T R^-1 q+eta^T R^-1 eta+q^T R^-1 H d+2 e^T J b-2 eta^T R^-1 H b+b^T J b+(Hb)^T R^-1(Hb)',
      'exact_rational_joseph_reduction_residual_zero':joseph_reduction==0,
      'exact_rational_reset_reduction_residual_zero':reset_reduction==0,
      'exact_rational_total_identity_residual_zero':original-reduced==0,
      'condition_number_scalarization_used':False,'independent_K_box_used':False,'independent_Sinverse_box_used':False,
      'posterior_precision_box_used':False,'packet_count_multiplier_used':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_shipping_P_H_R_K_cell_required','same_cell_cross_identity_consumed','S_inverse_eliminated_from_reduced_nonlinear_reset_ledger',
              'posterior_precision_eliminated_from_reduced_nonlinear_reset_ledger','prior_precision_is_same_moving_shipping_metric',
              'K_survives_only_through_same_cell_d_equals_Kq_graph','exact_rational_joseph_reduction_residual_zero',
              'exact_rational_reset_reduction_residual_zero','exact_rational_total_identity_residual_zero'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('condition_number_scalarization_used','independent_K_box_used','independent_Sinverse_box_used','posterior_precision_box_used','packet_count_multiplier_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
