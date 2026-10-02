# ruff: noqa: F401, F811
"""Blockwise accelerometer innovation Gram preserving the common rotation symbol.

For block-diagonal covariance at goLive:
 H=[-[f]x, 0..., R, I]
 HPH' = p_theta [f]x[f]x' + p_aw R R' + p_ba I
        = p_theta (||f||^2 I - f f') + (p_aw+p_ba) I.
Thus ALL R dependence cancels exactly before enclosure.
"""
from __future__ import annotations
import math,numpy as np
from .interval_riccati_21 import IMat
from .golive_release_seed import certificate as seed

def exact_golive_gram(f):
 f=np.asarray(f,float);c=seed();d=c["covariance_diagonal_upper"]
 pt=d[0];paw=d[15];pba=d[18]
 return pt*((f@f)*np.eye(3)-np.outer(f,f))+(paw+pba)*np.eye(3)

def force_cone_gram_bounds(rho_lo,rho_hi,axis,half_angle):
 """Entrywise outer Gram about cone-axis midpoint, after exact R cancellation."""
 axis=np.asarray(axis,float);axis/=np.linalg.norm(axis)
 rc=.5*(rho_lo+rho_hi);f0=rc*axis;G0=exact_golive_gram(f0)
 pt=seed()["covariance_diagonal_upper"][0]
 # A(f)=||f||^2 I-ff'. Lipschitz:
 # ||A(f)-A(g)||2 <= 2 (||f||+||g||) ||f-g||.
 # ||f-f0|| <= dr + rho_hi*2 sin(alpha/2).
 dr=.5*(rho_hi-rho_lo);df=dr+rho_hi*2*math.sin(half_angle/2)
 delta=2*pt*(rho_hi+rc)*df
 return G0,np.full((3,3),math.nextafter(delta,math.inf)),{
  "common_rotation_cancelled":True,"force_delta_norm":df,
  "gram_spectral_radius":delta}

def innovation_inverse_from_gram(Gmid,Grad,R):
 """Narrow inverse certificate from already-combined block Gram."""
 Rm=np.asarray(R.mid,float);Rr=np.asarray(R.rad,float)
 S0=Rm+np.asarray(Gmid,float);E=np.asarray(Grad,float)+Rr
 inv=np.linalg.inv(S0);eta=np.linalg.norm(inv,np.inf)*np.linalg.norm(E,np.inf)
 return {"verified":bool(eta<1),"residual":float(eta),
  "inverse_norm_bound":float(np.linalg.norm(inv,np.inf)/(1-eta)) if eta<1 else None,
  "common_rotation_cancelled":True}
