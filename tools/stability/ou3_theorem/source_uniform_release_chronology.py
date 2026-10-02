"""Literal source-uniform release chronology seam.

This module connects the explicit goLive seed, causal adaptation trace and
shipping S scheduler without inventing a finite magnetic callback pattern.
Magnetic service remains aggregate, exactly as the controlling theorem premise.

It is intentionally fail-closed: a source-uniform leaf is complete only when
same-history geometry/residual/local-defect enclosures are supplied.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from .golive_release_seed import certificate as seed_certificate
from .numeric_release_bounds import bounds as release_bounds

@dataclass
class Scheduler:
 elapsed:float=0.
 period:float=.015
 def retarget(self,new_period:float):
  if not(new_period>0 and math.isfinite(new_period)):raise ValueError("period")
  if self.elapsed>=new_period:self.elapsed=math.nextafter(new_period,0.)
  self.period=new_period
 def advance(self,dt:float):
  total=self.elapsed+dt
  tol=16*1.1920928955078125e-7*max(1.,self.period) # shipping float epsilon
  if total+tol<self.period:
   self.elapsed=total;return False
  self.elapsed=math.fmod(total,self.period) if total>=self.period else 0.
  if not(0<=self.elapsed<self.period):self.elapsed=0.
  return True

def pseudo_period(tau,ratio=.013636363636363634,lo=.005,hi=.15):
 return min(hi,max(lo,ratio*tau))

def release_seed():
 s=seed_certificate();b=release_bounds()
 return {"mean_state":s["mean_state"],
         "covariance_diagonal_upper":s["covariance_diagonal_upper"],
         "cross_covariance_zero":s["cross_covariance_zero_at_handoff"],
         "BA_graph_entry_interval":b["BA_graph_entry_interval"],
         "release_horizon_s":b["release_horizon_s"],
         "source_uniform_seed":True}

def schedule_from_adaptation(adaptation_trace,dt=.005,elapsed0=0.):
 """Enumerate literal S due events from causally generated applied tau."""
 sch=Scheduler(elapsed0,.015);out=[]
 for k,a in enumerate(adaptation_trace):
  # A non-point tau interval can straddle a due/not-due scheduler boundary.
  # Do not pick one branch in a promoted source-uniform leaf.
  plo=pseudo_period(a["tau"].lo);phi=pseudo_period(a["tau"].hi)
  if abs(phi-plo)>1e-15:
   out.append({"sample":k,"kind":"S_schedule_branch","period_lo":plo,
               "period_hi":phi,"resolved":False})
   continue
  sch.retarget(plo)
  if sch.advance(dt):out.append({"sample":k,"kind":"S","period":plo,"resolved":True})
 return out

def required_source_uniform_inputs():
 return ["attitude_rotation_geometry","acc_measurement_vector",
         "acc_local_defect_vector","physical_S_vector","estimator_S_vector",
         "accepted_acc_gate","reset_dtheta_from_same_Kr",
         "structured_covariance_coefficients"]

def connect(history_payload,geometry_stream=None):
 seed=release_seed();adapt=history_payload["adaptation_trace"]
 sched=schedule_from_adaptation(adapt)
 unresolved=[x for x in sched if not x["resolved"]]
 magnetic=history_payload["aggregate_magnetic_service"]
 missing=[] if geometry_stream is not None else required_source_uniform_inputs()
 return {"seed":seed,"S_schedule":sched,
         "S_schedule_source_uniform_resolved":not unresolved,
         "aggregate_magnetic_service":magnetic,
         "magnetic_callback_pattern_enumerated":False,
         "missing_same_history_inputs":missing,
         "ready_for_full_AbcD_iteration":not unresolved and not missing,
         "reason":None if not unresolved and not missing else
          "split shared history until S due branches resolve; attach same-history geometry/residual stream"}
