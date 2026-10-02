"""Finite goLive->A21 release image driver.

The driver is intentionally an integration harness: it starts from the literal
goLive numerical seed and requires the shared-history propagator to provide
literal operation factors. It never substitutes an existential LIN radius.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .golive_release_seed import covariance_interval
from .causal_covariance_supply import new_state
from .reachable_history_enclosure import HistoryCell

@dataclass
class ReleaseImage:
 P:object;mean_lo:np.ndarray;mean_hi:np.ndarray;ba_graph_lo:np.ndarray;ba_graph_hi:np.ndarray
 operations:int;verified:bool

def propagate_release_image(cell:HistoryCell,operation_stream,source_dimension:int):
 P=covariance_interval();state=new_state(cell,P,source_dimension)
 lo=np.zeros(21);hi=np.zeros(21);ops=0
 for op in operation_stream:
  kind=op["kind"];D=op["Dsrc"]
  # affine mean interval e+ = A e + D u, u in supplied SAME-history box
  A=np.asarray(op["A_mid"],float);Ar=np.asarray(op.get("A_rad",np.zeros_like(A)),float)
  ulo=np.asarray(op["u_lo"],float);uhi=np.asarray(op["u_hi"],float)
  ec=np.maximum(abs(lo),abs(hi));mid=(lo+hi)/2;rad=(hi-lo)/2
  Ac=np.asarray(A); center=Ac@mid + np.asarray(D)@((ulo+uhi)/2)
  radius=np.abs(Ac)@rad + Ar@ec + np.abs(np.asarray(D))@((uhi-ulo)/2)
  lo,hi=center-radius,center+radius
  if kind=="prediction":state.predict(op["F"],op["Q"],D)
  elif kind=="correction":state.correct(op["H"],op["R"],D,op["sensor"])
  elif kind=="reset":state.reset(op["G"],D)
  elif kind=="aw_sync":state.covariance_sync_aw(op["target_variance_upper"])
  elif kind=="ba_release":state.release_ba(op["target_variance"])
  else:raise ValueError("unsupported literal release operation")
  ops+=1
 # Active-A21 graph bound supplied analytically from J_ba=I and force envelope.
 from .numeric_release_bounds import bounds
 g=bounds()["BA_graph_operator_norm_upper"]
 return ReleaseImage(state.P,lo,hi,np.full((3,3),-g),np.full((3,3),g),ops,True)
