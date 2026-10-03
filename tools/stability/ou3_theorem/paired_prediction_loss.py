"""Prediction loss identity for paired covariance/probe storage.

For P-=F P F'+Q and Phi-=F Phi with invertible F, set C=F P F'.
Then the exact loss in Phi'P^-1Phi is
  X' [P^-1 - F'(C+Q)^-1 F] X
= Y' [C^-1-(C+Q)^-1] Y, Y=F X
= Y' C^-1 Q^(1/2) (I+Q^(1/2) C^-1 Q^(1/2))^-1 Q^(1/2) C^-1 Y.
This is PSD and is the only loss from a prediction step.
"""
from __future__ import annotations
import numpy as np

def prediction(P,X,F,Q):
 C=F@P@F.T; return C+Q,F@X,C

def loss(P,X,F,Q):
 Pn,Xn,C=prediction(P,X,F,Q)
 return X.T@np.linalg.solve(P,X)-Xn.T@np.linalg.solve(Pn,Xn)

def direct_loss(P,X,F,Q):
 _,Xn,C=prediction(P,X,F,Q)
 return Xn.T@(np.linalg.inv(C)-np.linalg.inv(C+Q))@Xn

def relative_bound(P,F,Q):
 """Scalar loss fraction bound: S_after >= S_before/(1+eta),
 eta=lambda_max(C^-1/2 Q C^-1/2)."""
 C=F@P@F.T
 L=np.linalg.cholesky(C); Li=np.linalg.inv(L)
 eta=float(np.linalg.eigvalsh((Li@Q@Li.T+Li@Q.T@Li.T)/2).max())
 return eta,1/(1+max(eta,0.0))

def certificate():
 rng=np.random.default_rng(4);A=rng.normal(size=(5,5));P=A@A.T+.5*np.eye(5)
 F=np.eye(5)+.02*rng.normal(size=(5,5));B=rng.normal(size=(5,5));Q=.001*(B@B.T)
 X=rng.normal(size=(5,2));L=loss(P,X,F,Q);D=direct_loss(P,X,F,Q)
 eta,retain=relative_bound(P,F,Q)
 return {"qualification":"OU3_PAIRED_PREDICTION_LOSS_IDENTITY_V1",
         "identity_defect":float(np.linalg.norm(L-D)),
         "loss_min_eigenvalue":float(np.linalg.eigvalsh((L+L.T)/2).min()),
         "eta":eta,"retention_factor":retain,
         "identity":"prediction is the only paired-storage loss; exact loss uses C^-1-(C+Q)^-1",
         "shipping_one_second_bound_verified":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
