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
from .interval_shaped_supply import IJointQuadratic,verified_precision

def _mid(a:IMat): return np.asarray(a.mid,float)
def correction_factor(k:IMat,h:IMat)->IMat:return add(eye(N),scale(matmul(k,h),-1))

@dataclass
class CovarianceSupplyState:
 root:HistoryCell
 P:IMat
 joint:JointQuadratic
 prefix: list
 interval_joint: IJointQuadratic|None=None
 def __post_init__(self):
  if self.P.shape!=(N,N):raise ValueError("A21 covariance required")
  if self.joint.dependency_token!=self.root.prefix_token:raise ValueError("dependency token mismatch")
  if self.interval_joint is None:self.interval_joint=IJointQuadratic.zeros(N,self.joint.B.shape[1],self.root.prefix_token)

 def _Qinterval(self):
  return verified_precision(self.P)

 def _Qmid(self):
  p=(_mid(self.P)+_mid(self.P).T)/2
  try:return np.linalg.inv(p)
  except np.linalg.LinAlgError as e:raise ArithmeticError("covariance midpoint storage singular") from e

 def predict(self,F:IMat,Q:IMat,Dsrc:np.ndarray):
  Q0,c0=self._Qinterval();q0=self._Qmid(); self.P=predict_covariance(self.P,F,Q);Q1,c1=self._Qinterval();q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(F),Dsrc)
  from .rank_loss_interval_factor import exact
  Di=exact(Dsrc.tolist())
  self.interval_joint.add_operation(Q0,Q1,F,Di);self.prefix.append({"mid":self.joint.matrix().copy(),"precision":[c0,c1]})

 def correct(self,H:IMat,R:IMat,Dsrc:np.ndarray,sensor:str):
  Q0,c0=self._Qinterval();q0=self._Qmid();S=innovation_covariance(self.P,H,R);Sinv,cert=verified_inverse3_interval(S)
  K=matmul(matmul(self.P,transpose(H)),Sinv);A=correction_factor(K,H)
  self.P=joseph_covariance(self.P,K,H,R);Q1,c1=self._Qinterval();q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(A),Dsrc)
  from .rank_loss_interval_factor import exact
  self.interval_joint.add_operation(Q0,Q1,A,exact(Dsrc.tolist()));self.prefix.append({"mid":self.joint.matrix().copy(),"precision":[c0,c1]})
  return {"sensor":sensor,"K":K,"A":A,"innovation_inverse":cert}

 def covariance_sync_aw(self,target_variance_upper:float):
  self.P=aw_covariance_floor_event(self.P,target_variance_upper)

 def release_ba(self,target_variance:float):
  self.P=accel_bias_release_event(self.P,target_variance)

 def reset(self,G:IMat,Dsrc:np.ndarray):
  # Shipping reset transports covariance congruently; no artificial Q.
  Q0,c0=self._Qinterval();q0=self._Qmid();self.P=matmul(matmul(G,self.P),transpose(G));Q1,c1=self._Qinterval();q1=self._Qmid()
  self.joint.add_operation(q0,q1,_mid(G),Dsrc)
  from .rank_loss_interval_factor import exact
  self.interval_joint.add_operation(Q0,Q1,G,exact(Dsrc.tolist()));self.prefix.append({"mid":self.joint.matrix().copy(),"precision":[c0,c1]})

def new_state(root:HistoryCell,P:IMat,source_dimension:int):
 return CovarianceSupplyState(root,P,JointQuadratic.zeros(N,source_dimension,root.prefix_token),[])

def source_columns_from_interval(mid:np.ndarray,rad:np.ndarray,token_count:int):
 """Embed one dependency-preserving affine source factor supplied by history propagator.

 The caller owns symbol identity; this function never invents independent
 columns for generated tuner/gain uncertainty.
 """
 mid=np.asarray(mid,float);rad=np.asarray(rad,float)
 if mid.shape!=(N,token_count) or rad.shape!=mid.shape:raise ValueError("source factor shape")
 if np.any(rad<0):raise ValueError("nonnegative source radii")
 return mid,rad
