"""Covariance-metric quotient by a certified one-dimensional compatibility line.

Analytical linear algebra only. The caller must supply the SAME-history
shipping covariance/precision and literal compatibility line. This module does
not assert that the line is admitted for all time under MAGNETIC SERVICE.
"""
from __future__ import annotations
import numpy as np

def whitened_line(P,line,tol=2e-10):
    """P=L L', y=L^-1 x; check projector identities in this J-isometry.

    An absolute residual in Pi' J-J Pi depends on the units/scale of P and
    falsely rejected the shipping word.  The same tolerance is applied to
    dimensionless whitened identities, without relaxing the SPD requirement.
    """
    P=np.asarray(P,float); line=np.asarray(line,float).reshape(-1)
    if P.ndim!=2 or P.shape[0]!=P.shape[1] or P.shape[0]!=line.size:
        raise ValueError("dimension mismatch")
    if not np.all(np.isfinite(P)) or not np.all(np.isfinite(line)):
        raise ValueError("finite covariance and line required")
    if not np.allclose(P,P.T,rtol=0,atol=tol*np.linalg.norm(P,2)):
        raise ValueError("symmetric covariance required")
    L=np.linalg.cholesky((P+P.T)/2)
    y=np.linalg.solve(L,line)
    norm=np.linalg.norm(y)
    if not (norm>0 and np.isfinite(norm)):
        raise ValueError("nonzero finite line in SPD metric")
    y=y/norm
    Q=np.outer(y,y); N=np.eye(line.size)-Q
    checks={"idempotent_Q":bool(np.linalg.norm(Q@Q-Q)<=tol),
            "idempotent_perp":bool(np.linalg.norm(N@N-N)<=tol),
            "J_orthogonal":bool(np.linalg.norm(Q.T-Q)<=tol),
            "annihilates_line":bool(np.linalg.norm(N@y)<=tol)}
    if not all(checks.values()): raise ArithmeticError("quotient projector identities")
    return L,y,Q,N,checks

def quotient(P,line,tol=2e-10):
    L,y,Q,N,checks=whitened_line(P,line,tol)
    W=np.linalg.solve(L,np.eye(len(y)))
    return W.T@W,L@y,L@Q@W,L@N@W,checks

def decompose(P,line,error):
    L,y,_,N,checks=whitened_line(P,line)
    e=np.asarray(error,float).reshape(-1)
    z=np.linalg.solve(L,e); zp=N@z
    alpha=float(y@z)
    V=float(z@z); Vp=float(zp@zp)
    defect=abs(V-(alpha*alpha+Vp))
    scale=max(1.,abs(V),abs(alpha*alpha+Vp))
    if defect>5e-10*scale: raise ArithmeticError("storage decomposition")
    return {"alpha":alpha,"V":V,"V_perp":Vp,"decomposition_defect":defect,
            "projector_checks":checks}

def quotient_terminal_map(P0,line0,PN,lineN,M,b):
    """Return transverse map and gauge-to-transverse injection in orthonormal J coordinates."""
    L0,y0,_,_,_=whitened_line(P0,line0); LN,yN,_,_,_=whitened_line(PN,lineN)
    from .kernel_restricted_action import householder_complement
    U0=householder_complement(y0); UN=householder_complement(yN)
    M=np.asarray(M,float); b=np.asarray(b,float).reshape(-1)
    if M.shape!=L0.shape or b.shape!=y0.shape or not np.all(np.isfinite(M)) or not np.all(np.isfinite(b)):
        raise ValueError("finite terminal map with matching dimension required")
    T=np.linalg.solve(LN,M@L0)
    MQ=UN.T@T@U0
    CQ=UN.T@T@y0
    bQ=UN.T@np.linalg.solve(LN,b)
    return {"M_Q":MQ,"C_Q":CQ,"b_Q":bQ,
            "gauge_injection_norm":float(np.linalg.norm(CQ)),
            "quotient_dimension":len(y0)-1}
