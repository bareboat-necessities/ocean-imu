"""Non-promoting algebra for ALT gravity quotient / magnetic-service words.

This module selects the replacement formulation, NOT a contraction certificate.
It does not generate physical histories, replace shipping state, or qualify a
word from labels. The canonical six supplied coordinates remain 18:24. A yaw
quotient is not permission to supply or discard axial gyro-bias error.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence
import numpy as np
QUALIFICATION="OU3_ALT_GRAVITY_QUOTIENT_MAGNETIC_SERVICE_V1"; JOINT_DIM=24; MOTION_DIM=18; NEUTRAL_INDICES=tuple(range(18,24)); FORMULATION="gravity_quotient_plus_informative_magnetic_superwords"
def selection_status():
 return {"qualification":QUALIFICATION,"formulation":FORMULATION,"phase":"feasibility_diagnostic","old_per_word_absolute_heading_target_retired":True,"old_ungauged_counterexample_retained":True,"no_mag_absolute_heading_contraction_required":False,"yaw_only_quotient_closes_axial_bias":False,"motion_coordinates_moved_to_bounded_supply":False,"canonical_supplied_indices":list(NEUTRAL_INDICES),"post_gauge_attempted_call_gap_s":0.04,"accepted_informative_gap_s":None,"transported_heading_pair_information_floor":None,"accepted_service_implied_by_call_schedule":False,"universal_first_magnetic_acquisition_claimed":False,"quotient_equivariance_closed":False,"uniform_endpoint_storage_coercivity_closed":False,"source_uniform_informative_service_closed":False,"corrected_source_uniform_rho":None,"corrected_rho_measured":False,"interval_enclosure_authorized":False,"storage_search_allowed":False,"ALT_STARTUP_PASS":False,"ALT_LIVE_PASS":False,"ALT_END_TO_END_PASS":False,"next_obligation":"measure a specified compatible storage ratio on same-history informative magnetic superwords; retain ungauged centre dynamics, all H18/A21 edges, and actual endpoint memory; no interval work yet"}
def validate_selection(report):
 e=selection_status(); return [f"selected formulation field differs: {k}" for k,v in e.items() if report.get(k)!=v]
def _matrix(v,n):
 a=np.asarray(v,dtype=float)
 if a.ndim!=2 or not np.isfinite(a).all(): raise ValueError(f"{n} must be a finite matrix")
 return a
def _spd(v,n):
 a=_matrix(v,n)
 if a.shape[0]!=a.shape[1] or not a.shape[0]: raise ValueError(f"{n} must be a nonempty square matrix")
 if not np.allclose(a,a.T,rtol=0,atol=1e-12): raise ValueError(f"{n} must be symmetric")
 a=(a+a.T)/2
 try:return a,np.linalg.cholesky(a)
 except np.linalg.LinAlgError as e: raise ValueError(f"{n} must be positive definite") from e
def coordinate_injection(n,keep):
 idx=tuple(keep)
 if n<=0 or not idx or len(set(idx))!=len(idx) or any(not isinstance(i,(int,np.integer)) or i<0 or i>=n for i in idx): raise ValueError("invalid coordinate injection")
 return np.eye(n)[:,idx]
def joint24_motion_injection(): return coordinate_injection(JOINT_DIM,tuple(range(MOTION_DIM)))
def projected_storage_ratio(A,M_before,M_after,Z):
 A=_matrix(A,"A"); Z=_matrix(Z,"Z"); before,_=_spd(M_before,"M_before"); after,_=_spd(M_after,"M_after")
 if A.shape!=(after.shape[0],before.shape[0]) or Z.shape[0]!=before.shape[0] or Z.shape[1]==0: raise ValueError("A/storage/kernel dimensions do not match")
 den,L=_spd(Z.T@before@Z,"restricted M_before"); AZ=A@Z; num=AZ.T@after@AZ; left=np.linalg.solve(L,num); W=np.linalg.solve(L,left.T).T; W=(W+W.T)/2; vals,vecs=np.linalg.eigh(W); y=vecs[:,-1]; u=np.linalg.solve(L.T,y); x=Z@u; r=float(vals[-1]); residual=num@u-r*den@u
 return {"rho_point":r,"distance_to_one":1-r,"maximizing_direction":x.tolist(),"generalized_eigen_residual_norm":float(np.linalg.norm(residual)),"storage_before_min_eigenvalue":float(np.linalg.eigvalsh(before)[0]),"storage_after_min_eigenvalue":float(np.linalg.eigvalsh(after)[0]),"restricted_storage_condition_number":float(np.linalg.cond(den)),"point_below_one":r<1,"roundoff_outward_enclosed":False,"source_uniform_rho_certified":False}
def common_storage_family_ratio(maps,M,Z):
 maps=list(maps)
 if not maps: raise ValueError("a family must contain at least one map")
 rows=[projected_storage_ratio(A,M,M,Z) for A in maps]; w=max(range(len(rows)),key=lambda i:rows[i]["rho_point"])
 return {"metric_policy":"one_common_metric","word_count":len(rows),"worst_word_index":w,"worst_rho_point":rows[w]["rho_point"],"rows":rows,"finite_family_is_universal_source_cover":False}
def quotient_map(A,N_before,N_after,Q_before,Q_after,tolerance=1e-12):
 if not math.isfinite(tolerance) or tolerance<0: raise ValueError("nonnegative finite tolerance required")
 A=_matrix(A,"A"); N0,N1,Q0,Q1=[_matrix(v,n) for v,n in ((N_before,"N_before"),(N_after,"N_after"),(Q_before,"Q_before"),(Q_after,"Q_after"))]
 for N,Q,n in ((N0,Q0,A.shape[1]),(N1,Q1,A.shape[0])):
  frame=np.column_stack((N,Q))
  if N.shape[0]!=n or Q.shape[0]!=n or Q.shape[1]==0 or frame.shape!=(n,n) or not np.allclose(frame.T@frame,np.eye(n),rtol=0,atol=1e-11): raise ValueError("quotient bases must form a complete orthonormal frame")
 leak=Q1.T@A@N0; ln=float(np.linalg.norm(leak,2)) if leak.size else 0
 if ln>tolerance: raise ValueError(f"map does not descend to this quotient: leakage={ln}")
 return {"A_quotient":Q1.T@A@Q0,"fibre_to_quotient_leakage":ln,"point_invariance_only":True,"physical_gauge_equivariance_certified":False}
@dataclass(frozen=True)
class MagneticServiceContract:
 max_informative_gap_s:float; min_heading_response:float; min_pair_information:float
 def __post_init__(self):
  for n in ("max_informative_gap_s","min_heading_response","min_pair_information"):
   v=getattr(self,n)
   if not math.isfinite(v) or v<=0: raise ValueError(f"{n} must be finite and positive")
@dataclass(frozen=True)
class MagneticEvent:
 time_s:float; accepted:bool; gauged:bool; transported_whitened_heading_rows:object=None
def audit_magnetic_service(events,start_s,end_s,contract):
 if not(math.isfinite(start_s) and math.isfinite(end_s) and start_s<end_s): raise ValueError("finite increasing window endpoints required")
 times=[]; info=np.zeros((2,2)); last=start_s; accepted=0
 for e in events:
  t=e.time_s
  if not math.isfinite(t) or t<start_s or t>end_s or t<last: raise ValueError("events must be chronological and inside the window")
  if type(e.accepted) is not bool or type(e.gauged) is not bool: raise ValueError("accepted and gauged must be literal booleans")
  last=t
  if not(e.accepted and e.gauged): continue
  accepted+=1
  if e.transported_whitened_heading_rows is None: raise ValueError("accepted gauged event lacks transported information")
  rows=_matrix(e.transported_whitened_heading_rows,"heading rows")
  if rows.shape[1]!=2 or rows.shape[0]==0: raise ValueError("heading rows must have two columns and positive row count")
  if float(np.linalg.norm(rows[:,0]))<contract.min_heading_response: continue
  times.append(t); info+=rows.T@rows
 gaps=np.diff([start_s,*times,end_s]); worst=float(np.max(gaps)); eig=np.linalg.eigvalsh(info); gp=bool(times) and worst<=contract.max_informative_gap_s; ip=bool(eig[0]>=contract.min_pair_information)
 return {"calls":len(events),"accepted_gauged_calls":accepted,"informative_calls":len(times),"worst_informative_gap_s":worst,"information_gramian":info.tolist(),"information_min_eigenvalue":float(eig[0]),"accepted_gap_point_pass":gp,"pair_information_point_pass":ip,"finite_window_point_service_pass":bool(gp and ip),"runtime_acceptance_guaranteed":False,"infinite_recurrence_certified":False,"source_uniform_observability_certified":False}
@dataclass(frozen=True)
class WordPiece:
 A:object; history_id:str; before_id:str; after_id:str
def compose_carried_pieces(pieces):
 pieces=list(pieces)
 if not pieces: raise ValueError("at least one carried piece required")
 first=_matrix(pieces[0].A,"A")
 if first.shape[0]!=first.shape[1] or first.shape[0]==0: raise ValueError("pieces must be nonempty square maps")
 result=np.eye(first.shape[0]); previous=None; history=pieces[0].history_id
 for p in pieces:
  A=_matrix(p.A,"A")
  if A.shape!=first.shape or not history or p.history_id!=history or not p.before_id or not p.after_id: raise ValueError("piece dimensions/history/boundary identifiers differ")
  if previous is not None and p.before_id!=previous: raise ValueError("superword has a detached/reseeded boundary")
  result=A@result; previous=p.after_id
 return result
