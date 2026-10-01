"""Fail-closed interval factor certificates for the OU-III rank-loss proof.

This module propagates literal *factor* maps, not scalar covariance envelopes.
Every source column remains shared across observations.  Midpoint/radius
operations use analytic radii plus nextafter expansion.  Matrix inverses are
accepted only after an exact-binary64 residual/Neumann certificate.

The module deliberately does not invent shipping boxes.  Callers must supply
interval F/U/H/V/A matrices from the literal event chronology.  If a source
box is too coarse, certification returns a zero lower bound with a reason.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
import math

from .interval_riccati_21 import IMat, add, matmul, scale, transpose
from .interval_riccati import symmetric_interval_gershgorin


def _out(x: float) -> float:
    return math.nextafter(x, math.inf)


def exact(rows) -> IMat:
    rows=tuple(tuple(float(x) for x in r) for r in rows)
    return IMat(rows,tuple(tuple(0.0 for _ in r) for r in rows))


def zeros(n: int,m: int) -> IMat:
    if n<=0 or m<0:
        raise ValueError("positive row and nonnegative column counts required")
    return exact([[0.0]*m for _ in range(n)])


def eye(n: int) -> IMat:
    return exact([[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)])


def hstack(a: IMat,b: IMat) -> IMat:
    if a.shape[0]!=b.shape[0]: raise ValueError("row mismatch")
    return IMat(tuple(x+y for x,y in zip(a.mid,b.mid)),
                tuple(x+y for x,y in zip(a.rad,b.rad)))


def vstack(a: IMat,b: IMat) -> IMat:
    if a.shape[1]!=b.shape[1]: raise ValueError("column mismatch")
    return IMat(a.mid+b.mid,a.rad+b.rad)


def sub(a: IMat,b: IMat) -> IMat:
    return add(a,scale(b,-1.0))


def block_diag(blocks: tuple[IMat,...]) -> IMat:
    nr=sum(x.shape[0] for x in blocks); nc=sum(x.shape[1] for x in blocks)
    mm=[[0.0]*nc for _ in range(nr)]; rr=[[0.0]*nc for _ in range(nr)]
    i0=j0=0
    for x in blocks:
        n,m=x.shape
        for i in range(n):
            for j in range(m):
                mm[i0+i][j0+j]=x.mid[i][j]; rr[i0+i][j0+j]=x.rad[i][j]
        i0+=n;j0+=m
    return IMat(tuple(tuple(x) for x in mm),tuple(tuple(x) for x in rr))


def _mid_inverse(a):
    n=len(a)
    x=[list(map(float,row))+[1.0 if i==j else 0.0 for j in range(n)]
       for i,row in enumerate(a)]
    for k in range(n):
        p=max(range(k,n),key=lambda i:abs(x[i][k]))
        if x[p][k]==0.0: raise ArithmeticError("singular midpoint")
        x[k],x[p]=x[p],x[k]
        q=x[k][k]
        x[k]=[z/q for z in x[k]]
        for i in range(n):
            if i==k: continue
            q=x[i][k]
            x[i]=[u-q*v for u,v in zip(x[i],x[k])]
    return tuple(tuple(row[n:]) for row in x)


def verified_inverse(a: IMat) -> tuple[IMat,dict]:
    """Generic square inverse enclosure from an exact binary64 residual test."""
    n,m=a.shape
    if n!=m: raise ValueError("square matrix required")
    try:
        b=_mid_inverse(a.mid)
    except ArithmeticError as error:
        return zeros(n,n),{"verified":False,"reason":str(error)}
    if any(not math.isfinite(x) for row in b for x in row):
        return zeros(n,n),{"verified":False,"reason":"non-finite midpoint inverse"}
    # Exact residual of the binary64 midpoint/inverse pair.
    rmax=Fraction(0)
    for i in range(n):
        row=Fraction(0)
        for j in range(n):
            s=sum((Fraction.from_float(b[i][k])*Fraction.from_float(a.mid[k][j])
                   for k in range(n)),Fraction(0))
            rij=(Fraction(1) if i==j else Fraction(0))-s
            row+=abs(rij)
        rmax=max(rmax,row)
    bnorm=max(sum(abs(Fraction.from_float(x)) for x in row) for row in b)
    erad=max(sum(Fraction.from_float(max(0.0,x)) for x in row) for row in a.rad)
    q=rmax+bnorm*erad
    if q>=1:
        try:
            q_upper=float(q)
        except OverflowError:
            q_upper=math.inf
        return zeros(n,n),{"verified":False,"residual_ratio_upper":q_upper}
    delta=bnorm*q/(1-q)
    try:
        radius=_out(float(delta))
    except OverflowError:
        radius=math.inf
    if not math.isfinite(radius):
        return zeros(n,n),{"verified":False,"reason":"non-finite inverse enclosure"}
    rad=tuple(tuple(radius for _ in range(n)) for _ in range(n))
    return IMat(b,rad),{"verified":True,"residual_ratio_upper":float(q),
                        "inverse_error_inf_upper":float(delta)}


@dataclass(frozen=True)
class FactorState:
    """x = X root + B source for one literal linearized chronology."""
    X: IMat
    B: IMat

    def predict(self,F: IMat,U: IMat) -> "FactorState":
        if F.shape[1]!=self.X.shape[0] or U.shape[0]!=F.shape[0]:
            raise ValueError("prediction shape mismatch")
        return FactorState(matmul(F,self.X),hstack(matmul(F,self.B),U))

    def affine_mean(self,A: IMat,U: IMat|None=None) -> "FactorState":
        out=FactorState(matmul(A,self.X),matmul(A,self.B))
        return out if U is None else FactorState(out.X,hstack(out.B,U))

    def observe(self,H: IMat,V: IMat) -> tuple[IMat,IMat]:
        """Return y=O root + C source, with fresh V columns appended."""
        o=matmul(H,self.X)
        c=hstack(matmul(H,self.B),V)
        return o,c

    def correct(self,A: IMat,K: IMat,V: IMat) -> "FactorState":
        """Literal e+=A e- - K V n; sign is immaterial to covariance/action."""
        kv=scale(matmul(K,V),-1.0)
        return FactorState(matmul(A,self.X),hstack(matmul(A,self.B),kv))


def source_covariance(c: IMat) -> IMat:
    # An empty factor has zero covariance; IMat does not represent zero rows.
    if c.shape[1]==0:
        return zeros(c.shape[0],c.shape[0])
    return matmul(c,transpose(c))


def schur_scalar_information(v0: IMat,va: IMat,residual_cov: IMat) -> dict:
    """Verified scalar Schur information of va after eliminating columns v0."""
    if va.shape[1]!=1 or v0.shape[0]!=va.shape[0]:
        raise ValueError("va must be one column with matching rows")
    rinv,rcert=verified_inverse(residual_cov)
    if not rcert["verified"]:
        return {"verified":False,"lower":0.0,"reason":"residual covariance inverse"}
    gaa=matmul(matmul(transpose(va),rinv),va)
    if v0.shape[1]==0:
        gam=gaa
        gcert={"verified":True,"reason":"no nuisance columns"}
    else:
        g00=matmul(matmul(transpose(v0),rinv),v0)
        g0a=matmul(matmul(transpose(v0),rinv),va)
        g00i,gcert=verified_inverse(g00)
        if not gcert["verified"]:
            return {"verified":False,"lower":0.0,"reason":"nuisance Gram inverse"}
        loss=matmul(matmul(transpose(g0a),g00i),g0a)
        gam=sub(gaa,loss)
    lo=math.nextafter(gam.mid[0][0]-gam.rad[0][0],-math.inf)
    hi=math.nextafter(gam.mid[0][0]+gam.rad[0][0],math.inf)
    return {"verified":lo>0.0,"lower":max(0.0,lo),"upper":hi,
            "residual_inverse":rcert,"nuisance_inverse":gcert}


def residualized_gram(obs: IMat,nuis: IMat,residual_cov: IMat) -> tuple[IMat,dict]:
    """O'R^-1O shorted by nuisance columns, with verified inverses."""
    rinv,rcert=verified_inverse(residual_cov)
    if not rcert["verified"]:
        return zeros(obs.shape[1],obs.shape[1]),{"verified":False,"lower":0.0,
                                                "reason":"residual covariance inverse"}
    goo=matmul(matmul(transpose(obs),rinv),obs)
    if nuis.shape[1]==0:
        g=goo
    else:
        gnn=matmul(matmul(transpose(nuis),rinv),nuis)
        gno=matmul(matmul(transpose(nuis),rinv),obs)
        gnni,ncert=verified_inverse(gnn)
        if not ncert["verified"]:
            return zeros(obs.shape[1],obs.shape[1]),{"verified":False,"lower":0.0,
                                                    "reason":"nuisance Gram inverse"}
        g=sub(goo,matmul(matmul(transpose(gno),gnni),gno))
    g=scale(add(g,transpose(g)),.5)
    lo,hi=symmetric_interval_gershgorin(g.mid,g.rad)
    return g,{"verified":lo>0.0,"lower":max(0.0,lo),"upper":hi,
              "residual_inverse":rcert}


def generalized_ratio_lower(action_gram: IMat,terminal_gram: IMat) -> dict:
    """Conservative verified lower bound min x'Ax/x'Tx.

    Gershgorin is intentionally fail-closed.  Adaptive callers should split
    boxes or change coordinates when this lower bound contains zero.
    """
    alo,ahi=symmetric_interval_gershgorin(action_gram.mid,action_gram.rad)
    tlo,thi=symmetric_interval_gershgorin(terminal_gram.mid,terminal_gram.rad)
    if thi<=0:
        return {"verified":True,"lower":math.inf,"reason":"terminal-null"}
    lower=max(0.0,alo)/thi
    return {"verified":alo>0.0 and thi>0.0,"lower":lower,
            "action_eigen_lower":alo,"terminal_eigen_upper":thi}


def source_range_audit() -> dict:
    """Current theorem-range audit: fail closed until literal factor boxes exist."""
    return {
        "qualification":"OU3_RANK_LOSS_INTERVAL_FACTOR_V1",
        "literal_ranges":{
            "tau_s":[0.02,12.0],"dt_s":[0.004,0.006],
            "S_gap_max_s":0.15,"sigma_aw_min":0.05,
            "accel_std_min":0.05,"integral_std_min":0.075,
            "aw_covariance_upper":16.0,"R_S_variance_upper":10000.0,
            "magnetic_service_window_s":1.0,"magnetic_service_mu":1.0,
        },
        "gamma_S_lower":0.0,
        "gamma_M_lower":0.0,
        "beta_lower":0.0,
        "source_uniform_verified":False,
        "reason":"literal common-source interval factor boxes not yet supplied",
    }
