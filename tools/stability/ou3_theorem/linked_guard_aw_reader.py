"""Backward AW reader through literal release events."""
from __future__ import annotations
import numpy as np
from .guard_weighted_fast_abel import integrated_primitive_bound

def mmid(M): return np.asarray(M.mid,float)

def terminal_aw_reader(events,n_samples,ell=(0.,0.,1.)):
 lam=np.zeros(21);lam[15:18]=np.asarray(ell,float)
 coeff=np.zeros((n_samples,3))
 for e in reversed(events):
  kind=e.get("kind")
  if kind=="acc":
   k=int(e["sample"]);K=mmid(e["K"])
   coeff[k]+=lam@K
   H=mmid(e["H"])
   lam=lam@(np.eye(21)-K@H)
  elif kind=="S":
   K=mmid(e["K"]);H=mmid(e["H"])
   lam=lam@(np.eye(21)-K@H)
  elif kind=="prediction":
   lam=lam@mmid(e["F"])
  # reset adjoint is already represented in event G.
  if "G" in e: lam=lam@mmid(e["G"])
 return coeff

def linked_charge(events,guard_weights,dt,C,H,alpha,ell=(0.,0.,1.)):
 n=len(guard_weights);R=terminal_aw_reader(events,n,ell)
 wm=[.5*(w.lo+w.hi) for w in guard_weights]
 wr=[.5*(w.hi-w.lo) for w in guard_weights]
 axes=[integrated_primitive_bound(alpha,dt,C,0.,H,wm,R[:,j]) for j in range(3)]
 return {"axis_bounds":[z["bound"] for z in axes],
         "axis_charges":[z["charge"] for z in axes],
         "max_midpoint_bound":max(z["bound"] for z in axes),
         "reader":R,"weight_mid":wm,"weight_rad":wr,
         "promotion_ready":False,
         "remaining":"interval adjoint reader and weight-radius contribution"}
