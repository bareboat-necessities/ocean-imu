# ruff: noqa: F401, F811
"""Causal vibration-conditioned Racc enclosure from the same guard history."""
from __future__ import annotations
import math,numpy as np
from .causal_tuner_interval import I

def effective_sigma(base,scales,excess,gain=.75):
 b=np.asarray(base,float);s=np.asarray(scales,float)
 return np.sqrt((b*s)**2+(gain*float(excess))**2)

def covariance_from_sigma(sigma):
 s=np.asarray(sigma,float);return np.diag(s*s)

def covariance_interval_from_excess(excess:I,base=.2,scales=(1.,1.,1.),gain=.75):
 s=np.asarray(scales,float);lo=(base*s)**2+(gain*excess.lo)**2
 hi=(base*s)**2+(gain*excess.hi)**2
 mid=.5*(lo+hi);rad=np.nextafter(.5*(hi-lo),np.inf)
 return {"mid":np.diag(mid),"rad":np.diag(rad),"linked":True,
         "same_guard_history":True}

def qualified_covariance_ceiling(std_max=.3010398644698074):
 return np.eye(3)*(std_max*std_max)

def source_uniform(history_excess=None,base=.2,scales=(1.,1.,1.)):
 if history_excess is None:
  return {"R":qualified_covariance_ceiling(),"linked":False,
          "qualification_ceiling_used":True,
          "reason":"same-history AccelVibrationGuard state not propagated"}
 if isinstance(history_excess,I):
  return covariance_interval_from_excess(history_excess,base,scales)
 lo,hi=history_excess
 return covariance_interval_from_excess(I(lo,hi),base,scales)
