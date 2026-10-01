"""Interval magnetic-service rows with explicit E_hb/nuisance split.

MAGNETIC SERVICE protects the two E_hb columns of the ACTUAL transported and
whitened magnetic rows.  This module never treats that restricted floor as a
floor after arbitrary nuisance projection.  It builds the joint block
[G_hb,G_n] so the Schur loss caused by non-E_hb root columns is explicit.
"""
from __future__ import annotations
import math
from .interval_riccati_21 import IMat,matmul,transpose,add,scale
from .rank_loss_interval_factor import verified_inverse,zeros
from .interval_riccati import symmetric_interval_gershgorin


def _cols(a: IMat,indices):
    return IMat(tuple(tuple(row[j] for j in indices) for row in a.mid),
                tuple(tuple(row[j] for j in indices) for row in a.rad))


def _hstack_events(rows):
    # Vertical stack, despite historical helper naming.
    if not rows: raise ValueError("at least one magnetic event")
    n=rows[0].shape[1]
    if any(x.shape[1]!=n for x in rows): raise ValueError("root width mismatch")
    return IMat(tuple(r for x in rows for r in x.mid),
                tuple(r for x in rows for r in x.rad))


def whitened_root_rows(H: IMat,S_actual: IMat,Phi: IMat) -> tuple[IMat,dict]:
    """Enclose S^{-1/2} H Phi through the information Gram.

    We avoid interval Cholesky.  Return raw B=H Phi plus S^{-1}; downstream
    Grams use B' S^{-1} B exactly.
    """
    sinv,cert=verified_inverse(S_actual)
    if not cert["verified"]:
        return zeros(H.shape[0],Phi.shape[1]),{"verified":False,"reason":"innovation inverse"}
    return matmul(H,Phi),{"verified":True,"S_inverse":sinv,"inverse_certificate":cert}


def service_schur_information(events,ehb_indices,root_nuisance_indices):
    """Schur-short non-E_hb root columns from a literal 1-s magnetic block.

    events: iterable of (H_interval,S_actual_interval,Phi_interval).
    The actual innovation covariance is block diagonal across measurement
    source columns only for the measurement-noise part; chronological shared
    process sources must be supplied in Phi/source transport before this
    function is promoted.  This function therefore accepts each actual S and
    forms the exact information sum used by MAGNETIC SERVICE, then exposes
    the algebraic root-column Schur loss.
    """
    ghh=None;ghn=None;gnn=None
    inverse=[]
    for H,S,Phi in events:
        B,cert=whitened_root_rows(H,S,Phi)
        if not cert["verified"]:
            return {"verified":False,"gamma_M_lower":0.0,"reason":cert["reason"]}
        W=cert["S_inverse"];Bh=_cols(B,ehb_indices);Bn=_cols(B,root_nuisance_indices)
        ah=matmul(matmul(transpose(Bh),W),Bh)
        an=matmul(matmul(transpose(Bn),W),Bn) if root_nuisance_indices else zeros(0,0)
        hn=matmul(matmul(transpose(Bh),W),Bn) if root_nuisance_indices else zeros(len(ehb_indices),0)
        ghh=ah if ghh is None else add(ghh,ah)
        if root_nuisance_indices:
            ghn=hn if ghn is None else add(ghn,hn)
            gnn=an if gnn is None else add(gnn,an)
        inverse.append(cert["inverse_certificate"])
    if ghh is None: raise ValueError("empty service block")
    raw_lo,raw_hi=symmetric_interval_gershgorin(ghh.mid,ghh.rad)
    if not root_nuisance_indices:
        return {"verified":raw_lo>0,"gamma_M_lower":max(0.0,raw_lo),
                "unshorted_lower":raw_lo,"unshorted_upper":raw_hi,
                "inverse_certificates":inverse}
    gnni,ncert=verified_inverse(gnn)
    if not ncert["verified"]:
        return {"verified":False,"gamma_M_lower":0.0,
                "unshorted_lower":raw_lo,"reason":"nuisance Gram singular/uncertain",
                "inverse_certificates":inverse}
    schur=add(ghh,scale(matmul(matmul(ghn,gnni),transpose(ghn)),-1.0))
    schur=scale(add(schur,transpose(schur)),.5)
    lo,hi=symmetric_interval_gershgorin(schur.mid,schur.rad)
    return {"verified":lo>0.0,"gamma_M_lower":max(0.0,lo),
            "gamma_M_upper":hi,"unshorted_lower":raw_lo,
            "unshorted_service_floor_preserved":raw_lo>=1.0,
            "nuisance_inverse":ncert,"inverse_certificates":inverse,
            "schur":schur}


def service_floor_only_contract():
    """Fail-closed contract when only MAGNETIC SERVICE mu=1 is known.

    If G=[G_h,G_n], G_h'G_h>=I does not imply
    G_h'(I-P_Gn)G_h>0.  This function makes that non-implication executable.
    """
    return {"verified":False,"gamma_M_lower":0.0,"mu_M":1.0,
            "needs_literal_non_Ehb_transport":True,
            "reason":"restricted service floor does not bound root-column Schur complement"}
