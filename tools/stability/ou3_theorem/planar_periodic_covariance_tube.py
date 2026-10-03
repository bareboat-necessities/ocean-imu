"""Planar periodic covariance tube certificate framework.

The tube is centered on a literal settled periodic covariance profile, not on
the generic A21 root floor. A phase-indexed profile P0[j] is certified by
Loewner radii r[j] when every literal one-step map sends the entire input ball
P0[j]+/-r[j]I into the next ball. Scheduler phase is an explicit branch.

This module contains the exact Riccati monotonicity tests used after the literal
phase profile and factor rows are exported. It fails closed without that export.
"""
from __future__ import annotations
import numpy as np

def loewner_margin(A,B):
    return float(np.linalg.eigvalsh(((A-B)+(A-B).T)/2).min())

def joseph(P,H,R):
    S=H@P@H.T+R
    K=P@H.T@np.linalg.inv(S)
    A=np.eye(P.shape[0])-K@H
    return (A@P@A.T+K@R@K.T + (A@P@A.T+K@R@K.T).T)/2

def predict(P,F,Q):
    X=F@P@F.T+Q
    return (X+X.T)/2

def correction_monotone_check(Plo,Phi,H,R):
    """Positive-noise Riccati correction is Loewner monotone."""
    return loewner_margin(joseph(Phi,H,R),joseph(Plo,H,R))

def tube_step(Pcenter,r,F,Q,corrections,next_center,next_r):
    n=Pcenter.shape[0];I=np.eye(n)
    lo=Pcenter-r*I;hi=Pcenter+r*I
    if np.linalg.eigvalsh(lo).min()<=0:return {"contained":False,"reason":"input lower face not SPD"}
    lo=predict(lo,F,Q);hi=predict(hi,F,Q)
    for H,R in corrections:
        lo=joseph(lo,H,R);hi=joseph(hi,H,R)
    inner=next_center-next_r*I;outer=next_center+next_r*I
    return {"contained":loewner_margin(lo,inner)>=0 and loewner_margin(outer,hi)>=0,
            "lower_margin":loewner_margin(lo,inner),
            "upper_margin":loewner_margin(outer,hi)}

def certificate():
    return {"qualification":"OU3_PLANAR_PERIODIC_COVARIANCE_TUBE_V1",
            "representation":"phase-indexed Loewner balls around literal settled planar covariance profile, separately 12/9 parity",
            "scheduler":"branch over every S=0 placement induced by elapsed in [0,T_S)",
            "map":"literal F/Q prediction followed by literal positive-noise correction rows and reset congruence",
            "monotonicity":"prediction and positive-noise Riccati correction are Loewner monotone; reset is congruence",
            "profile_export_verified":False,
            "all_scheduler_branches_contained":False,
            "periodic_lower_tube_verified":False,
            "prediction_retention_lower":None,
            "required_retention":1/7.024764605642482,
            "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
