"""Generalized prediction eta on anisotropic planar factors.

Given a certified P>=L L' at a root, prediction gives C=F P F'>=F L L' F'.
Hence eta=lambda_max(Q,C)<=lambda_max(Q,F L L' F'). This module evaluates the
matrix generalized eigenvalue without scalarizing L. It also composes a
sequence only when a certified factor is supplied at every prefix.
"""
from __future__ import annotations
import numpy as np
from .planar_anisotropic_factors import generalized_eta

def eta_from_root_factor(F,Q,L):
    FL=np.asarray(F,float)@np.asarray(L,float)
    return generalized_eta(Q,FL)

def retention_from_etas(etas):
    r=1.0
    for x in etas:
      if x<0 or not np.isfinite(x): raise ValueError("finite nonnegative eta")
      r/=1+x
    return r

def certificate():
    # Synthetic anisotropic regression: a tiny BG floor does not contaminate
    # attitude when Q respects the same coordinates.
    L=np.diag([4e-4,1e-5]);F=np.array([[1.,.005],[0.,1.]])
    Q=np.diag([9.1125e-9,5e-13])
    eta=eta_from_root_factor(F,Q,L)
    return {"qualification":"OU3_PLANAR_GENERALIZED_PREDICTION_ETA_V1",
            "synthetic_eta":eta,"synthetic_retention":retention_from_etas([eta]*200),
            "formula":"eta_k <= lambda_max(Q_k, F_k L_k L_k^T F_k^T)",
            "literal_prefix_factors_verified":False,
            "one_second_product_verified":False,
            "reason":"existing anisotropic joint-root factor is certified at recurring roots, not yet transported as a lower factor through every intervening correction",
            "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
