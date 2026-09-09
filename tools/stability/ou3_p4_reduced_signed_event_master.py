#!/usr/bin/env python3
"""Matrix form of the reduced exact same-cell Joseph+reset P4 ledger.

For one accepted Joseph/reset event the exact reduced identity is

  DeltaV = -q'R^-1 q + eta'R^-1 eta + q'R^-1 H d
           +2 e'J b -2 eta'R^-1 H b
           +b'J b +(H b)'R^-1(H b).

All coordinates are retained as linear maps of ONE augmented graph coordinate z:

  e=E z, q=Q z, eta=N z, d=D z, b=B z.

Here d=Kq is the FULL H18/A21 Joseph correction. Only E_theta d drives the
quaternion reset downstream; retaining full d is essential for the H d cross
term in the exact reduced identity.

This module returns the symmetric interval matrix L such that DeltaV=z'Lz.
No S^-1, posterior precision, independent K box, packet count, or reset norm
port appears. The caller must bind d=Kq and the chord/reset/projection graph by
same-cell IQCs before asking the outward LDLT to prove negativity.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_add,matrix_mul,matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_joseph_reset_reduced_signed_identity as ID

QUALIFICATION='OU3_P4_REDUCED_SIGNED_EVENT_AUGMENTED_MASTER_V2'


def _shape(A):return len(A),len(A[0]) if A else 0
def _scale(A,a):
    c=Interval.point(float(a));return [[c*x for x in row] for row in A]
def _bilinear(A,W,B):
    X=matrix_mul(matrix_mul(matrix_transpose(A),W),B)
    return matrix_symmetric_hull(_scale(matrix_add(X,matrix_transpose(X)),0.5))
def _quad(A,W):return matrix_symmetric_hull(matrix_mul(matrix_mul(matrix_transpose(A),W),A))
def _sum(*As):
    if not As:raise ValueError('at least one matrix required')
    out=As[0]
    for A in As[1:]:out=matrix_symmetric_hull(matrix_add(out,A))
    return out

def reduced_event_matrix(E,Q,N,D,B,J,H,Rinv):
    er,n=_shape(E);qr,nq=_shape(Q);nr,nn=_shape(N);dr,nd=_shape(D);br,nb=_shape(B)
    jr,jc=_shape(J);hr,hc=_shape(H);rr,rc=_shape(Rinv)
    if n==0 or nq!=n or nn!=n or nd!=n or nb!=n:raise ValueError('all maps must share augmented coordinate')
    if er not in (18,21) or br!=er or dr!=er or (jr,jc)!=(er,er):raise ValueError('E/D/B/J dimensions must be H18/A21')
    if (qr,nr)!=(3,3) or (hr,hc)!=(3,er) or (rr,rc)!=(3,3):raise ValueError('measurement dimensions must be three')
    HB=matrix_mul(H,B)
    return _sum(
      _scale(_quad(Q,Rinv),-1.0),
      _quad(N,Rinv),
      _bilinear(Q,matrix_mul(Rinv,H),D),
      _scale(_bilinear(E,J,B),2.0),
      _scale(_bilinear(N,matrix_mul(Rinv,H),B),-2.0),
      _quad(B,J),
      _quad(HB,Rinv))

def qvalue(A,z):
    s=Interval.point(0.0)
    for i,row in enumerate(A):
        for j,a in enumerate(row):s=s+Interval.point(float(z[i]))*a*Interval.point(float(z[j]))
    return s

def _selector(rows,n,offset):
    A=[[Interval.point(0.0) for _ in range(n)] for _ in range(rows)]
    for i in range(rows):A[i][offset+i]=Interval.point(1.0)
    return A

def build():
    ident=ID.build();bad=ID.validate(ident)
    if bad:raise RuntimeError('reduced identity prerequisite failed: '+repr(bad))
    # z=[e18,q3,eta3,d18,b18].
    n=60;E=_selector(18,n,0);Q=_selector(3,n,18);N=_selector(3,n,21);D=_selector(18,n,24);B=_selector(18,n,42)
    I=Interval.point
    J=[[I((2.0+i/20.0) if i==j else (0.01 if abs(i-j)==1 else 0.0)) for j in range(18)] for i in range(18)]
    H=[[I(0.0) for _ in range(18)] for _ in range(3)]
    for i in range(3):H[i][i]=I(1.0);H[i][15+i]=I(0.2)
    Rinv=[[I(4.0 if i==j else 0.0) for j in range(3)] for i in range(3)]
    L=reduced_event_matrix(E,Q,N,D,B,J,H,Rinv)
    z=[0.01*(i-20) for i in range(n)]
    e=z[:18];q=z[18:21];eta=z[21:24];d=z[24:42];b=z[42:60]
    def mv(A,x):return [sum(A[i][j].lo*x[j] for j in range(len(x))) for i in range(len(A))]
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    R=[[x.lo for x in row] for row in Rinv];JJ=[[x.lo for x in row] for row in J];HH=[[x.lo for x in row] for row in H]
    Hb=mv(HH,b);Hd=mv(HH,d)
    direct=(-dot(q,mv(R,q))+dot(eta,mv(R,eta))+dot(q,mv(R,Hd))
            +2*dot(e,mv(JJ,b))-2*dot(eta,mv(R,Hb))+dot(b,mv(JJ,b))+dot(Hb,mv(R,Hb)))
    val=qvalue(L,z);defect=max(abs(val.lo-direct),abs(val.hi-direct))
    finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in L for x in row)
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'reduced_same_cell_identity_consumed':True,'augmented_matrix_API_available':True,
      'full_Joseph_correction_d_retained':True,'reset_uses_attitude_selector_of_d_downstream':True,
      'S_inverse_used':False,'posterior_precision_used':False,'independent_K_box_used':False,
      'independent_reset_norm_port_used':False,'packet_count_multiplier_used':False,
      'same_cell_d_equals_Kq_required_downstream':True,'dense_q_eta_d_b_e_cross_terms_retained':True,
      'matrix_finite':finite,'algebra_smoke_quadratic_defect':defect,
      'source_uniform_event_LDLT_closed_here':False,'source_uniform_endpoint_LDLT_closed_here':False,
      'source_uniform_every_prefix_LDLT_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('reduced_same_cell_identity_consumed','augmented_matrix_API_available','full_Joseph_correction_d_retained',
              'reset_uses_attitude_selector_of_d_downstream','same_cell_d_equals_Kq_required_downstream',
              'dense_q_eta_d_b_e_cross_terms_retained','matrix_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('S_inverse_used','posterior_precision_used','independent_K_box_used','independent_reset_norm_port_used',
              'packet_count_multiplier_used','source_uniform_event_LDLT_closed_here','source_uniform_endpoint_LDLT_closed_here',
              'source_uniform_every_prefix_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if not float(d.get('algebra_smoke_quadratic_defect',math.inf))<1e-10:f.append('reduced event matrix algebra smoke failed')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'matrix':d['augmented_matrix_API_available'],'defect':d['algebra_smoke_quadratic_defect'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
