"""Literal covariance/gain propagation from one causal shipping history cell.

Generated P/K are internal outputs.  Every mean factor and source column is fed
immediately into one joint shaped-supply accumulator with the same dependency
token. Proof-side only.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .interval_riccati_21 import (IMat,N,predict_covariance,innovation_covariance,
                                  verified_inverse3_interval,joseph_covariance,
                                  matmul,transpose,add,scale,aw_covariance_floor_event,
                                  accel_bias_release_event)
from .rank_loss_interval_factor import eye
from .reachable_history_enclosure import HistoryCell
from .joint_shaped_supply import JointQuadratic

def _mid(a:IMat): return np.asarray(a.mid,float)
def correction_factor(k:IMat,h:IMat)->IMat:return add(eye(N),scale(matmul(k,h),-1))

@dataclass
class CovarianceSupplyState:
 root:HistoryCell
 P:IMat
 joint:JointQuadratic
 prefix: list
 def __post_init__(self):
  if self.P.shape!=(N,N):raise ValueError("A21 covariance required")
  if self.joint.dependency_token!=self.root.prefix_token:raise ValueError("dependency token mismatch")

 def _Qmid(self):
  p=(_mid(self.P)+_mid(self.P).T)/2
  try:return np.linalg.inv(p)
  except np.linalg.LinAlgError as e:raise ArithmeticError("covariance midpoint storage singular") from e

 def predict(self,F:IMat,Q:IMat,Dsrc:np.ndarray):
  q0=self._Qmid(); self.P=predict_covariance(self.P,F,Q);q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(F),Dsrc);self.prefix.append(self.joint.matrix().copy())

 def correct(self,H:IMat,R:IMat,Dsrc:np.ndarray,sensor:str):
  q0=self._Qmid();S=innovation_covariance(self.P,H,R);Sinv,cert=verified_inverse3_interval(S)
  K=matmul(matmul(self.P,transpose(H)),Sinv);A=correction_factor(K,H)
  self.P=joseph_covariance(self.P,K,H,R);q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(A),Dsrc);self.prefix.append(self.joint.matrix().copy())
  return {"sensor":sensor,"K":K,"A":A,"innovation_inverse":cert}

 def covariance_sync_aw(self,target_variance_upper:float):
  self.P=aw_covariance_floor_event(self.P,target_variance_upper)

 def release_ba(self,target_variance:float):
  self.P=accel_bias_release_event(self.P,target_variance)

 def reset(self,G:IMat,Dsrc:np.ndarray):
  # Shipping reset transports covariance congruently; no artificial Q.
  q0=self._Qmid();self.P=matmul(matmul(G,self.P),transpose(G));q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(G),Dsrc);self.prefix.append(self.joint.matrix().copy())

def new_state(root:HistoryCell,P:IMat,source_dimension:int):
 return CovarianceSupplyState(root,P,JointQuadratic.zeros(N,source_dimension,root.prefix_token),[])
