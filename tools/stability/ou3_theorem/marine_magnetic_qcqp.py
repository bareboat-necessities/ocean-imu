"""Nonlinear MARINE/MAGNETIC whole-box callbacks for linked QCQP.

TRUE means every trajectory represented by the history box satisfies the
premise stratum encoded by its witness metadata. FALSE means disjoint.
UNKNOWN means subdivision is required. No premise is inferred from a midpoint.
"""
from __future__ import annotations
from dataclasses import dataclass
import math,numpy as np

TRUE=True;FALSE=False;UNKNOWN=None

@dataclass(frozen=True)
class VectorBox:
 lo:np.ndarray;hi:np.ndarray
 def norm_upper(self):return float(np.linalg.norm(np.maximum(np.abs(self.lo),np.abs(self.hi))))
 def norm_lower(self):
  q=np.where((self.lo<=0)&(self.hi>=0),0,np.minimum(np.abs(self.lo),np.abs(self.hi)))
  return float(np.linalg.norm(q))

def norm_cap(v:VectorBox,cap):
 if v.norm_upper()<=cap:return TRUE
 if v.norm_lower()>cap:return FALSE
 return UNKNOWN

def difference_span_witness(a:VectorBox,b:VectorBox,minimum):
 # Certifies ||a-b|| >= minimum for every realization in the two endpoint boxes.
 d=VectorBox(a.lo-b.hi,a.hi-b.lo)
 if d.norm_lower()>=minimum:return TRUE
 if d.norm_upper()<minimum:return FALSE
 return UNKNOWN

def gravity_direction_span_witness(g0:VectorBox,g1:VectorBox,theta):
 # Unit-vector boxes supplied by causal attitude/physical propagator.
 # For all realizations, dot <= cos(theta) certifies angle >= theta.
 maxdot=0.;mindot=0.
 for i in range(3):
  vals=(g0.lo[i]*g1.lo[i],g0.lo[i]*g1.hi[i],g0.hi[i]*g1.lo[i],g0.hi[i]*g1.hi[i])
  maxdot+=max(vals);mindot+=min(vals)
 c=math.cos(theta)
 if maxdot<=c:return TRUE
 if mindot>c:return FALSE
 return UNKNOWN

@dataclass(frozen=True)
class MagneticEventBox:
 time_lo:float;time_hi:float;field:VectorBox;residual:VectorBox;informative:bool|None

def magnetic_event_qualifies(e,field_min,field_max,resmax):
 nlo=e.field.norm_lower();nhi=e.field.norm_upper();rhi=e.residual.norm_upper()
 if e.informative is False or nhi<field_min or nlo>field_max or e.residual.norm_lower()>resmax:return FALSE
 if e.informative is True and nlo>=field_min and nhi<=field_max and rhi<=resmax:return TRUE
 return UNKNOWN

def recurring_magnetic_service(events,start,end,T,field_min,field_max,resmax):
 # Every complete T window needs an applied informative event. We certify by
 # covering [start,end-T] with event admissible-time intervals shifted by T.
 good=[];unknown=False
 for e in events:
  q=magnetic_event_qualifies(e,field_min,field_max,resmax)
  if q is TRUE:good.append((e.time_hi-T,e.time_lo))
  elif q is UNKNOWN:unknown=True
 target0=start;target1=end-T
 if target1<=target0:return TRUE
 good.sort();x=target0
 for a,b in good:
  if b<x:continue
  if a>x+1e-12:return UNKNOWN if unknown else FALSE
  x=max(x,b)
  if x>=target1:return TRUE
 return UNKNOWN if unknown else FALSE

@dataclass
class MarineMagneticWitness:
 velocity:list[VectorBox];position:list[VectorBox];acceleration:list[VectorBox];jerk:list[VectorBox]
 gravity_windows:list[list[tuple[VectorBox,VectorBox]]]
 displacement_windows:list[list[tuple[VectorBox,VectorBox]]]
 magnetic_events:list[MagneticEventBox];start:float;end:float

def callback(w,constants):
 m=constants["marine_motion"];mag=constants["magnetic_service"];states=[]
 for seq,cap in ((w.velocity,m["V_max_mps"]),(w.position,m["P_max_m"]),(w.acceleration,m["A_max_mps2"]),(w.jerk,m["J_max_mps3"])):
  for v in seq:states.append(norm_cap(v,cap))
 if any(x is FALSE for x in states):return FALSE
 if any(x is UNKNOWN for x in states):return UNKNOWN
 # Witness pairs are generated for every complete 30-s moving window by the
 # history-cell builder. Each pair must certify the required diameter/span.
 for pairs in w.gravity_windows:
  qs=[gravity_direction_span_witness(a,b,m["attitude_excitation"]["theta_E_rad"]) for a,b in pairs]
  if any(q is TRUE for q in qs): pass
  elif any(q is UNKNOWN for q in qs): return UNKNOWN
  else: return FALSE
 for pairs in w.displacement_windows:
  qs=[difference_span_witness(a,b,m["displacement_excitation"]["P_E_m"]) for a,b in pairs]
  if any(q is TRUE for q in qs): pass
  elif any(q is UNKNOWN for q in qs): return UNKNOWN
  else: return FALSE
 return recurring_magnetic_service(w.magnetic_events,w.start,w.end,mag["T_M_s"],
          mag["field_norm_min_uT"],mag["field_norm_max_uT"],mag["measurement_residual_norm_max_uT"])
