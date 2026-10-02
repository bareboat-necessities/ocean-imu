"""Aggregate MAGNETIC SERVICE operator for shaped action.

Consumes the theorem premise sum G_i'G_i >= mu I directly. No callback schedule
or individual service row enumeration is introduced.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .rank_loss_interval_factor import exact
from .interval_riccati_21 import add,IMat

@dataclass(frozen=True)
class AggregateMagneticService:
 window_s:float;mu:float;service_dimension:int=2
 def __post_init__(self):
  if not(self.window_s>0 and self.mu>0 and self.service_dimension==2):raise ValueError("MAGNETIC SERVICE premise")

def service_action_lower(service:AggregateMagneticService,embed:np.ndarray,weight_lower:float):
 """Return root-coordinate PSD action implied by aggregate service Gramian.

 If each accepted magnetic correction contributes at least weight_lower*G'G in
 the service coordinates, aggregate premise yields weight_lower*mu*E'E.
 The caller must certify weight_lower from the SAME covariance/gain history.
 """
 E=np.asarray(embed,float)
 if E.shape[0]!=2:raise ValueError("2-D service embedding")
 if not(weight_lower>0):raise ArithmeticError("positive same-history magnetic loss weight required")
 A=(weight_lower*service.mu)*(E.T@E)
 return exact(A.tolist()),{"aggregate_service_used":True,"event_schedule_enumerated":False,
                           "action_weight_lower":weight_lower,"mu":service.mu}

def augment_action(A:IMat,mag_action:IMat)->IMat:return add(A,mag_action)
