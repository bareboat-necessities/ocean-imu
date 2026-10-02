"""Single HistoryCell -> ConstructiveLeaf causal propagator.

This is the integration seam. It outer-interpolates theorem-input knots,
propagates the generated adaptation state sample-by-sample, and fails closed
where literal covariance/geometry inputs are not yet derivable from the cell.
"""
from __future__ import annotations
import json,math,numpy as np
from pathlib import Path
from .history_interpolation import sample_history,fast_primitive_increment_outer
from .causal_tuner_interval import (I,MahonyBox,initial_wave_period_box,BandBox,
 VarianceBox,ShippingTunerBox,ClosedCausalAdaptationBox)
from .history_witness_builder import TimedVectorBox,AttitudeGravityBox,build_history_witness
from .marine_magnetic_qcqp import VectorBox,MagneticEventBox
ROOT=Path(__file__).resolve().parents[3]
C=json.loads((ROOT/"tools/stability/ou3_theorem/constants.json").read_text())

def vb(xs):return VectorBox(np.array([x.lo for x in xs]),np.array([x.hi for x in xs]))

def initial_adaptation():
 z=I(0,0)
 return ClosedCausalAdaptationBox(
  MahonyBox((I(1,1),z,z,z),(z,z,z),True),initial_wave_period_box(),
  BandBox(z,z,z,z,z),VarianceBox(z,z,z,z),
  ShippingTunerBox(I(1.1,1.1),I(.1,.1),I(.5,.5)))

def propagate_history_cell(cell,horizon_s=60.,dt=.005):
 """Return integration payload or fail closed.

 A complete ConstructiveLeaf additionally needs literal initial covariance,
 acc/mag measurement geometry and BA graph enclosures from construction/release.
 Those are not theorem-input knot coordinates and MUST be supplied by the
 certified release-state propagator; absence is reported rather than invented.
 """
 n=int(round(horizon_s/dt));times=[i*dt for i in range(n+1)]
 samples=sample_history(cell,times,C);adapt=initial_adaptation();trace=[]
 vel=[];pos=[];acc=[];jerk=[];grav=[];mag=[]
 prev_a=None
 for i,row in enumerate(samples):
  t=row["t"];slow_a=tuple(row[f"slow_accel_{a}"] for a in range(3));slow_g=tuple(row[f"slow_gyro_{a}"] for a in range(3))
  # Physical acceleration is constrained by v derivative but sparse knots do
  # not determine it. Use theorem envelope; causal p/v consistency remains a
  # side constraint and prevents promotion if unresolved.
  Amax=C["marine_motion"]["A_max_mps2"];avec=(I(-Amax,Amax),)*3
  fast_a=tuple(fast_primitive_increment_outer(cell,"accel",a,max(0,t-dt),t,C)* (1/dt) for a in range(3))
  fast_g=tuple(fast_primitive_increment_outer(cell,"gyro",a,max(0,t-dt),t,C)* (1/dt) for a in range(3))
  delivered_a=tuple(avec[a]+slow_a[a]+fast_a[a] for a in range(3));delivered_g=tuple(slow_g[a]+fast_g[a] for a in range(3))
  try:state=adapt.step(I(dt,dt),delivered_g,delivered_a,I(.12,.12))
  except ArithmeticError as e:
   raise ArithmeticError("causal adaptation enclosure unresolved at t=%.6f: %s"%(t,e))
  trace.append(state)
  vv=vb(tuple(row[f"physical_v_{a}"] for a in range(3)));pp=vb(tuple(row[f"physical_p_{a}"] for a in range(3)));aa=vb(avec)
  vel.append(TimedVectorBox(t,t,vv));pos.append(TimedVectorBox(t,t,pp));acc.append(TimedVectorBox(t,t,aa))
  if prev_a is None:jv=VectorBox(np.full(3,-C["marine_motion"]["J_max_mps3"]),np.full(3,C["marine_motion"]["J_max_mps3"]))
  else:jv=VectorBox(np.full(3,-C["marine_motion"]["J_max_mps3"]),np.full(3,C["marine_motion"]["J_max_mps3"]))
  jerk.append(TimedVectorBox(t,t,jv));prev_a=aa
  gv=vb(tuple(row[f"gravity_dir_{a}"] for a in range(3)));grav.append(AttitudeGravityBox(t,t,gv))
  # Sparse magnetic knots cannot establish accepted/informative callbacks.
  mv=vb(tuple(row[f"mag_field_{a}"] for a in range(3)))
  mag.append(MagneticEventBox(t,t,mv,VectorBox(np.full(3,-C["magnetic_service"]["measurement_residual_norm_max_uT"]),np.full(3,C["magnetic_service"]["measurement_residual_norm_max_uT"])),None))
 # Do not invent release covariance/BA graph. Return the causal trace so the
 # release-state enclosure can be connected explicitly.
 return {"root":cell,"samples":samples,"adaptation_trace":trace,
         "physical":{"velocity":vel,"position":pos,"acceleration":acc,"jerk":jerk,"gravity":grav,"magnetic":mag},
         "complete_constructive_leaf":False,
         "reason":"certified release covariance/BA graph and accepted magnetic-event strata not yet connected"}
