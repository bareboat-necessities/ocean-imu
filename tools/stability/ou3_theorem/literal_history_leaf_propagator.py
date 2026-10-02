"""Single HistoryCell -> ConstructiveLeaf causal propagator.

This is the integration seam. It outer-interpolates theorem-input knots,
propagates the generated adaptation state sample-by-sample, and fails closed
where literal covariance/geometry inputs are not yet derivable from the cell.
"""
from __future__ import annotations
import json,math,numpy as np
from pathlib import Path
from .history_interpolation import sample_history,fast_primitive_increment_outer,physical_acceleration_history
from .causal_tuner_interval import (I,MahonyBox,initial_wave_period_box,BandBox,
 VarianceBox,ShippingTunerBox,ClosedCausalAdaptationBox)
from .history_witness_builder import TimedVectorBox,AttitudeGravityBox,build_history_witness
from .marine_magnetic_qcqp import VectorBox,MagneticEventBox
from .source_uniform_release_chronology import connect as connect_release_chronology
from .causal_vibration_guard_interval import initial_guard
from .causal_racc import covariance_interval_from_excess
from .guard_fast_state_abel import fast_lp_minus_raw_bound
from .enclosure_failure import EnclosureFailure
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
 samples=sample_history(cell,times,C);physical_acc=physical_acceleration_history(cell,times,C);adapt=initial_adaptation();guard=initial_guard();trace=[];fast_acc=[];racc=[];conditioned_acc=[];guard_weights=[]
 vel=[];pos=[];acc=[];jerk=[];grav=[];mag=[]
 prev_a=None
 for i,row in enumerate(samples):
  t=row["t"];slow_a=tuple(row[f"slow_accel_{a}"] for a in range(3));slow_g=tuple(row[f"slow_gyro_{a}"] for a in range(3))
  # Same-history acceleration derived from physical velocity knots + jerk.
  avec=physical_acc[i]
  fast_a=tuple(fast_primitive_increment_outer(cell,"accel",a,max(0,t-dt),t,C)* (1.0/dt) for a in range(3))
  fast_g=tuple(fast_primitive_increment_outer(cell,"gyro",a,max(0,t-dt),t,C)* (1.0/dt) for a in range(3))
  fast_acc.append(fast_a)
  # Body specific-force components: gravity direction is carried explicitly.
  # Physical translational acceleration has only a world-frame norm contract at
  # this seam, so its body components are conservatively rotated into [-A,A].
  # The coupled norm floor below retains the reverse-triangle information lost
  # by that component hull.
  Amax=C["marine_motion"]["A_max_mps2"];g=9.80665
  gv_comp=tuple(row[f"gravity_dir_{a}"] for a in range(3))
  body_a=(I(-Amax,Amax),)*3
  delivered_a=tuple(g*gv_comp[a]+body_a[a]+slow_a[a]+fast_a[a] for a in range(3));delivered_g=tuple(slow_g[a]+fast_g[a] for a in range(3))
  acc_norm_floor=g-Amax-C["imu_bias"]["B_a_s_mps2"]-C["imu_bias"]["B_a_f_mps2"]
  raw_step=(C["marine_motion"]["J_max_mps3"]+g*C["marine_motion"]["Omega_max_rad_s"]+Amax*C["marine_motion"]["Omega_max_rad_s"]+C["imu_bias"]["D_a_s_mps3"])*dt
  alpha_guard=math.exp(-2*math.pi*guard.cutoff_hz*dt)
  fast_state=fast_lp_minus_raw_bound(alpha_guard,guard.poles,dt,C["imu_bias"]["B_a_f_mps2"],C["imu_bias"]["fast_accel_accumulation_cap_mps"],C["imu_bias"]["fast_accel_horizon_s"],t)
  gs=guard.step(delivered_a,dt,acc_norm_floor,raw_step,fast_state);guard_weights.append(gs["weight"]);conditioned_acc.append(gs["conditioned"]);racc.append(covariance_interval_from_excess(gs["excess"]))
  conditioned_norm_floor=gs["conditioned_norm_lower"]
  try:state=adapt.step(I(dt,dt),delivered_g,gs["conditioned"],I(.12,.12),acc_norm_lower=conditioned_norm_floor)
  except (ArithmeticError,ValueError,OverflowError) as e:
   sens={name:iv.width for name,iv in cell.coordinates if name.startswith(("gravity_dir_","slow_accel_","fast_accel_primitive_","physical_v_","physical_p_"))}
   raise EnclosureFailure("causal adaptation",str(e),t,sens) from e
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
 from .aggregate_magnetic_service import AggregateMagneticService
 mag_service=AggregateMagneticService(C["magnetic_service"]["T_M_s"],C["magnetic_service"]["mu_M"])
 payload={"root":cell,"samples":samples,"adaptation_trace":trace,"physical_acceleration":physical_acc,"fast_accel":fast_acc,"conditioned_accel":conditioned_acc,"Racc_interval":racc,"guard_weight":guard_weights,
         "physical":{"velocity":vel,"position":pos,"acceleration":acc,"jerk":jerk,"gravity":grav},
         "aggregate_magnetic_service":mag_service,
         "complete_constructive_leaf":False}
 release=connect_release_chronology(payload)
 payload["release_chronology"]=release
 payload["reason"]="same-history attitude/measurement geometry and residual/local-defect stream not yet attached; explicit goLive seed and aggregate magnetic service are available"
 return payload
