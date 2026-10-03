"""Covariance-metric quotient by a certified one-dimensional compatibility line.

Analytical linear algebra only.  The caller must supply the SAME-history
shipping covariance/precision and literal compatibility line.
"""
from __future__ import annotations
import numpy as np

def quotient(P,line,tol=2e-10):
    P=np.asarray(P,float); r=np.asarray(line,float).reshape(-1)
    if P.ndim!=2 or P.shape[0]!=P.shape[1] or P.shape[0]!=r.size:
        raise ValueError("dimension mismatch")
    J=np.linalg.inv(P)
    q=float(r@J@r)
    if not(q>0 and np.isfinite(q)): raise ValueError("nonzero finite line in SPD metric")
    r=r/np.sqrt(q)
    PiQ=np.outer(r,r@J); Pperp=np.eye(r.size)-PiQ
    checks={"idempotent_Q":np.linalg.norm(PiQ@PiQ-PiQ)<=tol,
            "idempotent_perp":np.linalg.norm(Pperp@Pperp-Pperp)<=tol,
            "J_orthogonal":np.linalg.norm(PiQ.T@J-J@PiQ)<=tol,
            "annihilates_line":np.linalg.norm(Pperp@r)<=tol}
    if not all(checks.values()): raise ArithmeticError("quotient projector identities")
    return J,r,PiQ,Pperp,checks

def decompose(P,line,error):
    J,r,PiQ,Pperp,checks=quotient(P,line)
    e=np.asarray(error,float).reshape(-1)
    alpha=float(r@J@e); ep=Pperp@e
    V=float(e@J@e); Vp=float(ep@J@ep)
    defect=abs(V-(alpha*alpha+Vp))
    scale=max(1.,abs(V),abs(alpha*alpha+Vp))
    if defect>5e-10*scale: raise ArithmeticError("storage decomposition")
    return {"alpha":alpha,"V":V,"V_perp":Vp,"decomposition_defect":defect,
            "projector_checks":checks}

def quotient_terminal_map(P0,line0,PN,lineN,M,b):
    """Return transverse map and gauge-to-transverse injection in orthonormal J coordinates."""
    J0,r0,_,_,_=quotient(P0,line0); JN,rN,_,_,_=quotient(PN,lineN)
    # Cholesky J=L L^T; Euclidean complements in y=L^T e coordinates.
    L0=np.linalg.cholesky(J0); LN=np.linalg.cholesky(JN)
    y0=L0.T@r0; yN=LN.T@rN
    from .kernel_restricted_action import householder_complement
    U0=householder_complement(y0); UN=householder_complement(yN)
    E0=np.linalg.solve(L0.T,U0); EN=np.linalg.solve(LN.T,UN)
    M=np.asarray(M,float); b=np.asarray(b,float).reshape(-1)
    MQ=UN.T@LN.T@M@E0
    CQ=UN.T@LN.T@M@r0
    bQ=UN.T@LN.T@b
    return {"M_Q":MQ,"C_Q":CQ,"b_Q":bQ,
            "gauge_injection_norm":float(np.linalg.norm(CQ)),
            "quotient_dimension":len(r0)-1}
