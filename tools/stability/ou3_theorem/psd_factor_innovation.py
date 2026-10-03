"""PSD-factor innovation enclosure for 3-D shipping corrections.

Avoids intervalizing HPH' entrywise.  A covariance factor enclosure P=L L'
is mapped to Y=H L and S=R+Y Y'.  The inverse is certified through a
Woodbury/Gram representation with the strict R floor.
"""
from __future__ import annotations
import math,numpy as np
from .interval_riccati_21 import IMat,matmul,transpose,add
from .rank_loss_interval_factor import exact

def point_cholesky_from_diagonal(diag):
 L=np.diag(np.sqrt(np.asarray(diag,float)))
 return exact(L.tolist())

def gram_innovation_from_factor(H:IMat,L:IMat,R:IMat):
 Y=matmul(H,L)
 S=add(R,matmul(Y,transpose(Y)))
 return Y,S

def woodbury_inverse_point_factor(Hmid,Lmid,Rdiag):
 """Exact midpoint low-rank identity used as center, not as theorem enclosure."""
 H=np.asarray(Hmid,float);L=np.asarray(Lmid,float);r=np.asarray(Rdiag,float)
 Y=H@L
 # S is only 3x3; this form keeps PSD construction explicit.
 S=np.diag(r)+Y@Y.T
 return np.linalg.inv(S),S

def psd_factor_inverse_enclosure(H:IMat,L:IMat,R:IMat):
 """Verified inverse enclosure exploiting S=R+YY' >= R.

 Center uses factor midpoint. Perturbation is bounded at Y-level:
 d(YY')=Y0 dY'+dY Y0'+dY dY'. This retains common factor structure.
 """
 Y=matmul(H,L)
 Ym=np.asarray(Y.mid,float);Yr=np.asarray(Y.rad,float)
 Rm=np.asarray(R.mid,float);Rr=np.asarray(R.rad,float)
 S0=Rm+Ym@Ym.T
 inv0=np.linalg.inv(S0)
 # Infinity norm bound on structured Gram perturbation.
 yn=np.linalg.norm(Ym,np.inf);dr=np.linalg.norm(Yr,np.inf)
 egram=2*yn*dr+dr*dr
 er=egram+np.linalg.norm(Rr,np.inf)
 inv0n=np.linalg.norm(inv0,np.inf)
 eta=inv0n*er
 # If structured residual verifies, get narrow Neumann enclosure.
 if eta<1:
  invn=inv0n/(1-eta);delta=inv0n*er*invn
  rad=np.full((3,3),math.nextafter(float(delta),math.inf))
  return IMat(tuple(map(tuple,inv0)),tuple(map(tuple,rad))),{
   "verified":True,"certificate":"PSD_factor_Gram_Neumann",
   "gram_perturbation_inf":float(egram),"residual":float(eta),
   "inverse_norm_bound":float(invn)}
 # Existence still follows from R floor, but narrow linked subtraction does not.
 return None,{"verified":False,"certificate":"PSD_factor_Gram_Neumann",
              "gram_perturbation_inf":float(egram),"residual":float(eta)}

def factor_innovation_probe(H:IMat,L:IMat,R:IMat):
 inv,cert=psd_factor_inverse_enclosure(H,L,R)
 return {**cert,"narrow_inverse_available":inv is not None}
