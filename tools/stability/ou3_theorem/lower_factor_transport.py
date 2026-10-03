"""Transport a covariance lower factor through a positive-noise correction.

If P>=L L' and measurement information J=H'R^-1H<=Jbar, inverse order gives
P+ = (P^-1+J)^-1 >= ( (L L')^-1 + Jbar )^-1.
This supplies a matrix-valued posterior lower factor without scalarization.
"""
from __future__ import annotations
import numpy as np

def corrected_lower(L,Jbar):
 L=np.asarray(L,float);J=np.asarray(Jbar,float)
 B=L@L.T
 C=np.linalg.inv(np.linalg.inv(B)+J)
 C=(C+C.T)/2
 return np.linalg.cholesky(C)

def check(P,L,H,R,Jbar):
 S=H@P@H.T+R;K=P@H.T@np.linalg.inv(S);A=np.eye(P.shape[0])-K@H
 Pp=A@P@A.T+K@R@K.T;Lp=corrected_lower(L,Jbar)
 return float(np.linalg.eigvalsh((Pp-Lp@Lp.T+Pp.T-Lp@Lp.T)/2).min())

def certificate():
 rng=np.random.default_rng(12);L=np.diag([.4,.03,.2]);Z=rng.normal(size=(3,3));P=L@L.T+Z@Z.T
 H=rng.normal(size=(2,3));R=.2*np.eye(2);J=H.T@np.linalg.solve(R,H)
 Lp=corrected_lower(L,J)
 return {"qualification":"OU3_MATRIX_LOWER_FACTOR_CORRECTION_TRANSPORT_V1",
         "posterior_factor":Lp.tolist(),"synthetic_loewner_margin":check(P,L,H,R,J),
         "identity":"P>=LL', J<=Jbar => Pplus >= ((LL')^-1+Jbar)^-1",
         "literal_information_ceilings_verified":False,
         "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
