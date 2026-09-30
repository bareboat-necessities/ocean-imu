"""Literal OU-III interval operation boxes used by magnetic event export.

These constructors mirror shipping formulas.  They are proof-side only and
never alter estimator behavior.  Every constructor is fail-closed; adaptive
callers must split state/covariance boxes when a verified inverse or gate
stratum cannot be established.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from .interval_riccati_21 import IMat,add,matmul,scale,transpose
from .rank_loss_interval_factor import eye,zeros,verified_inverse
from .magnetic_literal_box_export import PredictionBox,CorrectionBox,ResetBox


def _skew(v: IMat)->IMat:
    if v.shape!=(3,1): raise ValueError("3-vector interval required")
    m=[r[0] for r in v.mid];q=[r[0] for r in v.rad]
    mm=((0.,-m[2],m[1]),(m[2],0.,-m[0]),(-m[1],m[0],0.))
    rr=((0.,q[2],q[1]),(q[2],0.,q[0]),(q[1],q[0],0.))
    return IMat(mm,rr)


def magnetic_H_box(v2hat: IMat,nx:int=21)->IMat:
    """Literal shipping H_m: attitude block -[v2hat]x, all other columns zero."""
    J=scale(_skew(v2hat),-1.0)
    mm=[[0.0]*nx for _ in range(3)];rr=[[0.0]*nx for _ in range(3)]
    for i in range(3):
        for j in range(3):
            mm[i][j]=J.mid[i][j];rr[i][j]=J.rad[i][j]
    return IMat(tuple(tuple(x) for x in mm),tuple(tuple(x) for x in rr))


def magnetic_update_box(P:IMat,v2hat:IMat,Rmag:IMat,*,applied:bool,
                        t_lo=0.,t_hi=0.) -> tuple[CorrectionBox,dict]:
    """Construct H,S_actual,K,A=I-KH exactly as shipping mag-only update."""
    n=P.shape[0]
    if P.shape[1]!=n or Rmag.shape!=(3,3): raise ValueError("shape")
    H=magnetic_H_box(v2hat,n)
    S=add(matmul(matmul(H,P),transpose(H)),Rmag)
    if not applied:
        return CorrectionBox(eye(n),"mag",False,None,None,t_lo,t_hi),{
            "verified":True,"applied":False}
    Sinv,cert=verified_inverse(S)
    if not cert["verified"]:
        return CorrectionBox(eye(n),"mag",True,H,S,t_lo,t_hi),{
            "verified":False,"reason":"S_actual inverse","inverse":cert}
    K=matmul(matmul(P,transpose(H)),Sinv)
    A=add(eye(n),scale(matmul(K,H),-1.0))
    return CorrectionBox(A,"mag",True,H,S,t_lo,t_hi),{
        "verified":True,"applied":True,"K":K,"H":H,"S_actual":S,
        "inverse":cert}


def reset_box(dtheta:IMat,nx:int=21)->ResetBox:
    """Literal left-error reset mean/covariance transport G=I+.5[dtheta]x."""
    G=eye(nx);J=scale(_skew(dtheta),.5)
    mm=[list(r) for r in G.mid];rr=[list(r) for r in G.rad]
    for i in range(3):
        for j in range(3):
            mm[i][j]+=J.mid[i][j];rr[i][j]+=J.rad[i][j]
    return ResetBox(IMat(tuple(tuple(x) for x in mm),tuple(tuple(x) for x in rr)))


def prediction_box(F:IMat)->PredictionBox:
    return PredictionBox(F)


def covariance_after_correction(P:IMat,box:CorrectionBox,K:IMat,R:IMat)->IMat:
    """Joseph covariance image for operation-box propagation."""
    A=box.A
    return add(matmul(matmul(A,P),transpose(A)),
               matmul(matmul(K,R),transpose(K)))


@dataclass(frozen=True)
class MagneticStateBox:
    P: IMat
    Phi: IMat


def split_interval_matrix(a:IMat,i:int,j:int):
    """Bisect one interval entry; used by adaptive covariance/state subdivision."""
    r=a.rad[i][j]
    if r<=0: return (a,)
    m=a.mid[i][j];lo=m-r;hi=m+r;mid=.5*(lo+hi)
    outs=[]
    for xlo,xhi in ((lo,mid),(mid,hi)):
        mm=[list(x) for x in a.mid];rr=[list(x) for x in a.rad]
        mm[i][j]=.5*(xlo+xhi);rr[i][j]=math.nextafter(.5*(xhi-xlo),math.inf)
        outs.append(IMat(tuple(tuple(x) for x in mm),tuple(tuple(x) for x in rr)))
    return tuple(outs)


def widest_entry(a:IMat):
    best=(0,0,-1.0)
    for i,row in enumerate(a.rad):
        for j,r in enumerate(row):
            if r>best[2]:best=(i,j,r)
    return best


def adaptive_magnetic_update(P:IMat,v2hat:IMat,Rmag:IMat,*,applied=True,
                             max_depth=12,t_lo=0.,t_hi=0.):
    """Split the widest P/v2hat box until S inverse/gain verifies."""
    leaves=[]
    def rec(Pb,vb,depth):
        box,cert=magnetic_update_box(Pb,vb,Rmag,applied=applied,t_lo=t_lo,t_hi=t_hi)
        if cert["verified"]:
            leaves.append({"verified":True,"box":box,"cert":cert});return
        if depth>=max_depth:
            leaves.append({"verified":False,"box":box,"cert":cert});return
        ip,jp,rp=widest_entry(Pb);iv,jv,rv=widest_entry(vb)
        # Normalize vector width by a conservative 100 uT field scale.
        if rv/100.0 > rp/max(1.0,abs(Pb.mid[ip][jp])):
            for x in split_interval_matrix(vb,iv,jv):rec(Pb,x,depth+1)
        else:
            for x in split_interval_matrix(Pb,ip,jp):rec(x,vb,depth+1)
    rec(P,v2hat,0)
    return {"verified":all(x["verified"] for x in leaves),"leaves":leaves,
            "leaf_count":len(leaves)}
