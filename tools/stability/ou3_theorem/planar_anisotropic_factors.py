"""Assemble anisotropic planar parity lower factors from certified components.

This deliberately keeps matrix factors instead of scalarizing. The currently
available full-state Loewner certificate is the 16-s joint root comparison:
P >= diag(q_AG I6, A_LIN^-1 tensor I3, q_BA I3)/2.
It is permuted exactly into planar E/O coordinates and generalized Q/L
eigenvalues can then be evaluated without replacing L by min eig(L) I.
"""
from __future__ import annotations
import numpy as np
from .lin_matrix_certificate import action_matrix
from .root_covariance_certificate import process_floors
from .planar_parity import EVEN,ODD

def full_root_lower():
    a=np.array([[float(x) for x in row] for row in action_matrix()])
    lin=np.linalg.inv(a)/2.0
    ag,ba=process_floors(); ag=float(ag)/2.0;ba=float(ba)/2.0
    P=np.zeros((21,21));P[:6,:6]=ag*np.eye(6)
    # LIN state order is per-axis v,p,S,aw, while full order groups xyz blocks.
    for axis in range(3):
      idx=[6+axis,9+axis,12+axis,15+axis]
      P[np.ix_(idx,idx)]=lin
    P[18:,18:]=ba*np.eye(3)
    return P

def parity_factors():
    P=full_root_lower()
    out={}
    for name,idx in (("even",EVEN),("odd",ODD)):
      B=P[np.ix_(idx,idx)]
      L=np.linalg.cholesky(B)
      out[name]={"indices":list(idx),"lower_matrix":B.tolist(),"cholesky":L.tolist(),
                 "lambda_min":float(np.linalg.eigvalsh(B).min()),
                 "lambda_max":float(np.linalg.eigvalsh(B).max())}
    return out

def generalized_eta(Q,L):
    Q=np.asarray(Q,float);L=np.asarray(L,float)
    Li=np.linalg.inv(L);W=Li@Q@Li.T
    return float(np.linalg.eigvalsh((W+W.T)/2).max())

def certificate():
    f=parity_factors()
    return {"qualification":"OU3_PLANAR_ANISOTROPIC_PARITY_FACTORS_V1",
            "source":"existing exact 16-s joint-root Loewner certificate, permuted without scalarization",
            "even":f["even"],"odd":f["odd"],
            "generalized_eta_available":True,
            "one_second_literal_prediction_product_verified":False,
            "limitation":"factor is certified only at regular post-prediction roots after the 16-s path; corrections/predictions require propagation of the matrix lower factor to every one-second prefix",
            "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
