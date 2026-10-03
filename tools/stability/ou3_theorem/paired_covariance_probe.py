"""Exact paired covariance/probe identity for one Kalman correction.

For P+ = P-PH'(HPH'+R)^-1 HP and Phi+=(I-KH)Phi, the storage
Phi' P^-1 Phi is invariant across the correction. This is the correct joint
quantity; ordering P alone or Phi alone is generally invalid.
"""
from __future__ import annotations
import numpy as np

def correction(P,Phi,H,R):
 P=np.asarray(P,float);Phi=np.asarray(Phi,float);H=np.asarray(H,float);R=np.asarray(R,float)
 S=H@P@H.T+R;K=P@H.T@np.linalg.inv(S);A=np.eye(P.shape[0])-K@H
 Pp=A@P@A.T+K@R@K.T
 return Pp,A@Phi

def storage(P,Phi): return Phi.T@np.linalg.solve(P,Phi)

def identity_defect(P,Phi,H,R):
 Pp,Xp=correction(P,Phi,H,R)
 return float(np.linalg.norm(storage(Pp,Xp)-storage(P,Phi)))

def measurement_loss(Phi,H,S):
 Y=H@Phi
 return Y.T@np.linalg.solve(S,Y)

def prior_to_posterior_drop(P,Phi,H,R):
 S=H@P@H.T+R; Pp,Xp=correction(P,Phi,H,R)
 # The magnetic service summand is not the drop in Phi'P^-1Phi; that storage is invariant.
 return {"storage_defect":float(np.linalg.norm(storage(Pp,Xp)-storage(P,Phi))),
         "service":measurement_loss(Phi,H,S)}

def certificate():
 rng=np.random.default_rng(7);A=rng.normal(size=(6,6));P=A@A.T+.3*np.eye(6)
 X=rng.normal(size=(6,2));H=rng.normal(size=(3,6));R=.4*np.eye(3)
 d=identity_defect(P,X,H,R)
 return {"qualification":"OU3_PAIRED_COVARIANCE_PROBE_IDENTITY_V1",
         "identity_defect":d,
         "identity":"Phi_plus^T P_plus^-1 Phi_plus = Phi_minus^T P_minus^-1 Phi_minus for exact positive-noise Kalman/Joseph correction",
         "consequence":"nonmag accel/S corrections do not destroy paired covariance-metric probe storage; they redistribute it between covariance and homogeneous columns",
         "warning":"this invariant alone does not lower-bound the next magnetic innovation action Hm Phi because prediction/process noise and measurement geometry intervene",
         "shipping_service_proved":False}

if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
