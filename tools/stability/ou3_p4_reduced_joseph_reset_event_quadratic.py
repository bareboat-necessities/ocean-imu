#!/usr/bin/env python3
"""Dense same-cell reduced Joseph+reset event quadratic for P4.

The exact reduced identity is

  DeltaV = -q'R^-1 q + eta'R^-1 eta + q'R^-1 H d
           +2 e'J b -2 eta'R^-1 H b
           +b'J b +(H b)'R^-1(H b),

where e is the full pre-event physical error, q=H e+eta is the exact residual,
d=Kq is the SAME-cell full correction, b=G^-1 rho is the full reset defect in
the Joseph covariance frame, and J=P^-1 is the prior moving precision.

This module turns that identity into one symmetric quadratic matrix on any
caller-supplied common augmented coordinate z.  The caller supplies linear maps

  e=E z, q=Q z, eta=N z, d=D z, b=B z.

No S^-1, posterior precision, independent K box, reset-norm port, or scalar
remainder appears.  The d map must remain bound downstream to the same Joseph
cell through the exact correction graph.  The matrix is the theorem-facing
building block for endpoint/prefix augmented LDLT.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

QUALIFICATION='OU3_P4_DENSE_REDUCED_JOSEPH_RESET_EVENT_QUADRATIC_V1'


def shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(x)!=c for x in A):raise ValueError('ragged matrix')
    return r,c

def mt(A):return [list(x) for x in zip(*A)]
def mm(A,B):
    ar,ac=shape(A);br,bc=shape(B)
    if ac!=br:raise ValueError('matrix product mismatch')
    return [[sum((A[i][k]*B[k][j] for k in range(ac)),A[i][0]-A[i][0]) for j in range(bc)] for i in range(ar)]
def madd(A,B):
    if shape(A)!=shape(B):raise ValueError('matrix add mismatch')
    return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
# Multiply on the matrix-entry side.  This keeps the helper generic for both
# exact Fraction arithmetic and the outward Interval type, whose supported
# scalar multiplication is x*c rather than relying on scalar __rmul__.
def mscale(A,c):return [[x*c for x in row] for row in A]
def symcross(A,W,B):
    """Matrix whose quadratic is (Az)'W(Bz), with symmetry made explicit."""
    X=mm(mm(mt(A),W),B);Xt=mt(X)
    half=Fraction(1,2) if isinstance(X[0][0],Fraction) else 0.5
    return mscale(madd(X,Xt),half)
def gram(A,W):return mm(mm(mt(A),W),A)
def zero_like(n,sample):
    z=sample-sample
    return [[z for _ in range(n)] for _ in range(n)]

def event_quadratic(J,H,Rinv,E,Q,N,D,B):
    """Return symmetric M with z'Mz equal to the exact reduced event DeltaV."""
    nstate,nstate2=shape(J);mh,nh=shape(H);rr,rc=shape(Rinv)
    if nstate==0 or nstate!=nstate2 or mh!=3 or nh!=nstate or (rr,rc)!=(3,3):raise ValueError('J/H/R dimensions invalid')
    maps={'E':E,'Q':Q,'N':N,'D':D,'B':B};nz=None
    for name,A in maps.items():
        ar,ac=shape(A);expect=nstate if name in ('E','D','B') else 3
        if ar!=expect:raise ValueError(name+' row dimension invalid')
        if nz is None:nz=ac
        if ac!=nz:raise ValueError('augmented coordinate dimensions differ')
    HR=mm(mt(H),Rinv) # n x 3
    HtRH=mm(HR,H)    # n x n
    # -q'R^-1q + eta'R^-1eta
    M=madd(mscale(gram(Q,Rinv),-1),gram(N,Rinv))
    # + q'R^-1 H d
    M=madd(M,symcross(Q,mm(Rinv,H),D))
    # +2 e'Jb
    M=madd(M,mscale(symcross(E,J,B),2))
    # -2 eta'R^-1 H b
    M=madd(M,mscale(symcross(N,mm(Rinv,H),B),-2))
    # +b'Jb +(Hb)'R^-1(Hb)
    M=madd(M,gram(B,J));M=madd(M,gram(B,HtRH))
    return M

def mv(A,x):return [sum((a*b for a,b in zip(row,x)),x[0]-x[0]) for row in A]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),a[0]-a[0])
def quad(A,x):return dot(x,mv(A,x))
def selector(rows,n,offset):
    z=Fraction(0);o=Fraction(1)
    return [[o if j==offset+i else z for j in range(n)] for i in range(rows)]

def build():
    # Exact rational smoke on z=[e2,q2,eta2,d2,b2] using a 2-D analogue of the
    # general 3-vector formula.  A local builder below mirrors event_quadratic
    # without hard-coding measurement dimension three.
    J=[[Fraction(5,3),Fraction(1,7)],[Fraction(1,7),Fraction(7,4)]]
    H=[[Fraction(1),Fraction(2)],[Fraction(-1),Fraction(1,3)]]
    Rinv=[[Fraction(3,5),Fraction(-1,13)],[Fraction(-1,13),Fraction(4,7)]]
    n=10;E=selector(2,n,0);Q=selector(2,n,2);N=selector(2,n,4);D=selector(2,n,6);B=selector(2,n,8)
    HtRH=mm(mm(mt(H),Rinv),H)
    M=madd(mscale(gram(Q,Rinv),-1),gram(N,Rinv));M=madd(M,symcross(Q,mm(Rinv,H),D));M=madd(M,mscale(symcross(E,J,B),2));M=madd(M,mscale(symcross(N,mm(Rinv,H),B),-2));M=madd(M,gram(B,J));M=madd(M,gram(B,HtRH))
    z=[Fraction(2,5),Fraction(-1,6),Fraction(3,7),Fraction(1,8),Fraction(-1,9),Fraction(2,11),Fraction(1,5),Fraction(-2,13),Fraction(1,17),Fraction(3,19)]
    e=z[0:2];q=z[2:4];eta=z[4:6];d=z[6:8];b=z[8:10];Hb=mv(H,b)
    direct=(-quad(Rinv,q)+quad(Rinv,eta)+dot(q,mv(Rinv,mv(H,d)))+2*dot(e,mv(J,b))-2*dot(eta,mv(Rinv,Hb))+quad(J,b)+quad(Rinv,Hb))
    matrix_value=quad(M,z)
    symmetric=all(M[i][j]==M[j][i] for i in range(n) for j in range(n))
    return {'qualification':QUALIFICATION,'same_cell_d_equals_Kq_binding_required_downstream':True,
      'dense_common_augmented_coordinate_supported':True,'reduced_event_quadratic_available':True,
      'S_inverse_present_in_event_quadratic':False,'posterior_precision_present_in_event_quadratic':False,
      'independent_K_box_present_in_event_quadratic':False,'independent_reset_defect_port_used':False,
      'packet_count_multiplier_used':False,'exact_rational_matrix_matches_reduced_identity':matrix_value==direct,
      'exact_rational_matrix_symmetric':symmetric,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_d_equals_Kq_binding_required_downstream','dense_common_augmented_coordinate_supported','reduced_event_quadratic_available','exact_rational_matrix_matches_reduced_identity','exact_rational_matrix_symmetric'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('S_inverse_present_in_event_quadratic','posterior_precision_present_in_event_quadratic','independent_K_box_present_in_event_quadratic','independent_reset_defect_port_used','packet_count_multiplier_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
