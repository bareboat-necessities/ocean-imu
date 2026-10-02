"""Causal vibration-conditioned Racc enclosure.

A full interval replica of AccelVibrationGuard is not yet available.  This
module therefore exposes the exact shipping algebra once the same-history
excess-RMS interval is supplied, and a rigorous qualified ceiling fallback.
The fallback is source-uniform but deliberately loses correlation with r_acc.
"""
from __future__ import annotations
import math,numpy as np

def effective_sigma(base,scales,excess,gain=.75):
 b=np.asarray(base,float);s=np.asarray(scales,float)
 return np.sqrt((b*s)**2+(gain*float(excess))**2)

def covariance_from_sigma(sigma):
 s=np.asarray(sigma,float);return np.diag(s*s)

def qualified_covariance_ceiling(std_max=.3010398644698074):
 return np.eye(3)*(std_max*std_max)

def source_uniform(history_excess=None,base=.2,scales=(1.,1.,1.)):
 if history_excess is None:
  return {"R":qualified_covariance_ceiling(),"linked":False,
          "qualification_ceiling_used":True,
          "reason":"same-history AccelVibrationGuard state not propagated"}
 lo,hi=history_excess
 if lo!=hi:raise ArithmeticError("vibration excess branch not point-resolved")
 sig=effective_sigma(np.full(3,base),np.asarray(scales),lo)
 return {"R":covariance_from_sigma(sig),"linked":True,
         "qualification_ceiling_used":False}
