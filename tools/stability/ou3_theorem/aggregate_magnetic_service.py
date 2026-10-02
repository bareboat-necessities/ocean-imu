"""Aggregate MAGNETIC SERVICE information action.

The theorem rows are already transported and whitened by the actual innovation
covariance. Therefore sum G'G >= mu I is itself the homogeneous information
action; no extra gain weight or callback enumeration is introduced.
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

def service_action_lower(service:AggregateMagneticService,embed:np.ndarray):
 E=np.asarray(embed,float)
 if E.shape[0]!=2:raise ValueError("2-D service embedding")
 A=service.mu*(E.T@E)
 return exact(A.tolist()),{"aggregate_service_used":True,"event_schedule_enumerated":False,
  "actual_innovation_whitened":True,"action_weight_lower":1.0,"mu":service.mu}

def augment_action(A:IMat,mag_action:IMat)->IMat:return add(A,mag_action)
