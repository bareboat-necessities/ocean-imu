"""Exact moving-information prediction supply inequality.

For actual true-minus-estimate prediction error

 e_plus=F e+w,  P_plus=F P F^T+Q,

with P,Q positive definite,

 e_plus^T P_plus^-1 e_plus <= e^T P^-1 e+w^T Q^-1 w.

No invertibility of F is required. Let z=[e;w], D=diag(P,Q), B=[F I].
Then P_plus=B D B^T and
B^T(B D B^T)^-1 B <= D^-1 because the congruent matrix is an orthogonal
projector. Thus source/model/roundoff mismatch is an explicit supply channel
with no Young-factor loss. COMPLETE-BRMM/BIAS ancestry must still bind the same
history w downstream.
"""
from __future__ import annotations
from fractions import Fraction as Fq
import sys
from pathlib import Path
STABILITY=Path(__file__).resolve().parents[1]
if str(STABILITY) not in sys.path:sys.path.insert(0,str(STABILITY))
import ou3_full_process_ucc as PROCESS
from . import generalized_joseph as G
QUALIFICATION='OU3_ALT_MOVING_INFORMATION_PREDICTION_SUPPLY_V1'

def blockdiag(A,B):
    n=len(A);m=len(B);z=[[Fq(0) for _ in range(n+m)] for _ in range(n+m)]
    for i in range(n):
        for j in range(n):z[i][j]=A[i][j]
    for i in range(m):
        for j in range(m):z[n+i][n+j]=B[i][j]
    return z

def hstack(A,B):return [list(a)+list(b) for a,b in zip(A,B)]

def regression():
    P=[[Fq(3),Fq(1,4)],[Fq(1,4),Fq(2)]];Q=[[Fq(2),Fq(1,5)],[Fq(1,5),Fq(5,2)]];Fm=[[Fq(1),Fq(2,3)],[Fq(0),Fq(4,5)]]
    e=[Fq(2,7),Fq(-1,5)];w=[Fq(1,11),Fq(-2,13)]
    Pp=G.add(G.mm(G.mm(Fm,P),G.mt(Fm)),Q);ep=G.vadd(G.mv(Fm,e),w)
    lhs=G.quad(G.inverse(Pp),ep);rhs=G.quad(G.inverse(P),e)+G.quad(G.inverse(Q),w)
    D=blockdiag(P,Q);B=hstack(Fm,G.eye(2));gap=G.sub(G.inverse(D),G.mm(G.mm(G.mt(B),G.inverse(Pp)),B))
    A=[row[:] for row in gap];piv=[]
    for k in range(4):
        p=A[k][k];piv.append(p)
        if p<0:break
        if p==0:continue
        for i in range(k+1,4):
            l=A[i][k]/p
            for j in range(k+1,4):A[i][j]-=l*A[k][j]
    return {'lhs':lhs,'rhs':rhs,'inequality_exact':lhs<=rhs,'gap_pivots':piv,'gap_all_pivots_nonnegative':all(x>=0 for x in piv)}

def build():
    proc=PROCESS.build();pf=PROCESS.validate(proc)
    if pf:raise RuntimeError('process UCC prerequisite failed: '+repr(pf))
    h=float(proc['modes']['H']['prediction_Q_lambda_min_lower']);a=float(proc['modes']['A']['prediction_Q_lambda_min_lower']);r=regression()
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','exact_block_metric_prediction_supply_theorem':True,'F_invertibility_required':False,'P_positive_definite_required':True,'Q_positive_definite_required':True,'H18_shipping_Q_lambda_min_lower':h,'A21_shipping_Q_lambda_min_lower':a,'H18_Q_strictly_positive':h>0,'A21_Q_strictly_positive':a>0,'supply':'w_prediction^T Q^-1 w_prediction','Young_factor_or_epsilon_used':False,'same_history_prediction_forcing_must_be_bound_downstream':True,'source_or_model_forcing_declared_zero':False,'roundoff_declared_zero':False,'exact_rational_regression':r,'prediction_supply_structure_closed':r['inequality_exact'] and r['gap_all_pivots_nonnegative'] and h>0 and a>0,'source_uniform_prediction_forcing_bound_closed_here':False,'finite_precision_bound_closed_here':False,'source_uniform_complete_word_dissipation_closed_here':False,'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('exact_block_metric_prediction_supply_theorem','P_positive_definite_required','Q_positive_definite_required','H18_Q_strictly_positive','A21_Q_strictly_positive','same_history_prediction_forcing_must_be_bound_downstream','prediction_supply_structure_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('F_invertibility_required','Young_factor_or_epsilon_used','source_or_model_forcing_declared_zero','roundoff_declared_zero','source_uniform_prediction_forcing_bound_closed_here','finite_precision_bound_closed_here','source_uniform_complete_word_dissipation_closed_here','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if not float(d.get('H18_shipping_Q_lambda_min_lower',0))>0 or not float(d.get('A21_shipping_Q_lambda_min_lower',0))>0:f.append('shipping Q lower invalid')
    return f
