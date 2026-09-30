"""Aggregate MAGNETIC SERVICE Schur/canonical-correlation certificate.

No event schedule and no full covariance box are used.  The module records
the exact implication of the restricted service Gram and the additional
relative process modulus needed for a positive nuisance-shorted bound.
"""
from __future__ import annotations
import math

def service_only_counterexample():
    # G=[[I,c],[c',1]], c=(1,0): PSD rank two, service block I, Schur diag(0,1).
    return {"G_hh":[[1.,0.],[0.,1.]],"G_hn":[[1.],[0.]],"G_nn":[[1.]],
            "canonical_correlation":1.0,
            "schur":[[0.,0.],[0.,1.]],
            "service_mu":1.0,"gamma_residualized":0.0,
            "verified":True}

def joint_process_modulus(q_rel_lower:float):
    """If Q_n >= q_rel G_nn on the correlated nuisance range."""
    if q_rel_lower<=0 or not math.isfinite(q_rel_lower):
        return {"verified":False,"gamma_joint_lower":0.0,
                "reason":"positive relative nuisance process modulus required"}
    gamma=q_rel_lower/(1.0+q_rel_lower)
    return {"verified":True,"q_rel_lower":q_rel_lower,
            "eta_joint_upper":math.sqrt(1.0/(1.0+q_rel_lower)),
            "gamma_joint_lower":gamma}

def current_contract():
    return {"qualification":"OU3_MAG_SERVICE_AGGREGATE_SCHUR_V1",
            "service_only":service_only_counterexample(),
            "standalone_gamma_M_lower":0.0,
            "joint_q_rel_verified":False,
            "joint_gamma_lower":0.0,
            "next_obligation":"prove Q_n >= q_rel G_nn on Range(G_nh)"}


def relative_process_limit(M:float,q:float):
    if not (M>0 and q>0): raise ValueError("M,q positive")
    return {"M":M,"q":q,"service_mu":1.0,"q_rel_best":q/M,
            "joint_gamma":q/(M+q)}

def relative_process_route_status():
    return {"verified_negative":True,"source_uniform_q_rel_lower":0.0,
            "source_uniform_joint_gamma_lower":0.0,
            "reason":"service gives no upper or relative bound on correlated nuisance information",
            "next_obligation":"complete two-word Gram minimization with restricted service block retained"}
