"""Temporal same-history source domain for constructive OU-III supply.

Constraints attach to existing source symbols; no independent component maxima.
"""
from __future__ import annotations
from dataclasses import dataclass
import json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
C=json.loads((ROOT/"tools/stability/ou3_theorem/constants.json").read_text())

@dataclass(frozen=True)
class Symbol:
 name:str; kind:str; axis:int; time_s:float

@dataclass
class TemporalSourceDomain:
 symbols:tuple[Symbol,...]
 A:list[np.ndarray]; b:list[float]; notes:list[str]
 def add_norm_component_box(self,idxs,cap):
  # Necessary component inequalities; norm coupling is retained separately.
  for i in idxs:
   r=np.zeros(len(self.symbols));r[i]=1;self.A += [r,-r];self.b += [cap,cap]
 def add_difference(self,i,j,cap):
  r=np.zeros(len(self.symbols));r[i]=1;r[j]=-1;self.A += [r,-r];self.b += [cap,cap]
 def add_window_sum(self,idxs,weights,cap):
  r=np.zeros(len(self.symbols))
  for i,w in zip(idxs,weights):r[i]+=w
  self.A += [r,-r];self.b += [cap,cap]
 def matrices(self):
  return np.asarray(self.A,float),np.asarray(self.b,float)

def build_temporal_domain(symbols:tuple[Symbol,...],horizon_s:float):
 d=TemporalSourceDomain(symbols,[],[],[])
 imu=C["imu_bias"]; marine=C["marine_motion"]
 by=lambda kind,axis:[i for i,s in enumerate(symbols) if s.kind==kind and s.axis==axis]
 # SLOW: one carried path, amplitude plus pairwise rate/difference caps.
 for kind,B,D in (("slow_accel",imu["B_a_s_mps2"],imu["D_a_s_mps3"]),
                  ("slow_gyro",imu["B_g_s_rad_s"],imu["D_g_s_rad_s2"])):
  for a in range(3):
   ids=by(kind,a);d.add_norm_component_box(ids,B)
   for ii,i in enumerate(ids):
    for j in ids[ii+1:]:
     dt=abs(symbols[j].time_s-symbols[i].time_s);d.add_difference(i,j,min(2*B,D*dt))
 # FAST: instantaneous envelope and EVERY placed discrete window <= min(B*T,C).
 for kind,B,H,cap in (("fast_accel",imu["B_a_f_mps2"],imu["fast_accel_horizon_s"],imu["fast_accel_accumulation_cap_mps"]),
                       ("fast_gyro",imu["B_g_f_rad_s"],imu["fast_gyro_horizon_s"],imu["fast_gyro_accumulation_cap_rad"])):
  for a in range(3):
   ids=by(kind,a);d.add_norm_component_box(ids,B)
   for p in range(len(ids)):
    for q in range(p,len(ids)):
     T=max(0.,symbols[ids[q]].time_s-symbols[ids[p]].time_s)
     if T<=H+1e-12 and q>p:
      ws=[];sel=[]
      for k in range(p,q):
       dt=symbols[ids[k+1]].time_s-symbols[ids[k]].time_s
       sel.append(ids[k]);ws.append(dt)
      d.add_window_sum(sel,ws,min(B*T,cap))
 # Physical boundary variables: preserve bounded v,p and all-time acceleration.
 for kind,cap in (("physical_v",marine["V_max_mps"]),("physical_p",marine["P_max_m"]),("physical_a",marine["A_max_mps2"])):
  for a in range(3):d.add_norm_component_box(by(kind,a),cap)
 d.notes += ["Euclidean norm coupling and MARINE span are nonlinear side constraints, not replaced by component boxes.",
             "FAST windows are carried across proof boundaries; builder expects one persistent symbol timeline.",
             "assembled SLOW+FAST qualification remains conditional."]
 return d

def nonlinear_side_constraints(symbols):
 imu=C["imu_bias"];marine=C["marine_motion"];mag=C["magnetic_service"]
 return {"slow_norm_caps":{"accel":imu["B_a_s_mps2"],"gyro":imu["B_g_s_rad_s"]},
  "fast_norm_caps":{"accel":imu["B_a_f_mps2"],"gyro":imu["B_g_f_rad_s"]},
  "marine":{"velocity_norm":marine["V_max_mps"],"position_norm":marine["P_max_m"],
            "acceleration_norm":marine["A_max_mps2"],"jerk_norm":marine["J_max_mps3"],
            "attitude_span":{"T":marine["attitude_excitation"]["T_E_s"],"theta":marine["attitude_excitation"]["theta_E_rad"]},
            "position_span":{"T":marine["displacement_excitation"]["T_P_s"],"P":marine["displacement_excitation"]["P_E_m"]}},
  "magnetic_service":{"T":mag["T_M_s"],"mu":mag["mu_M"],"residual_norm":mag["measurement_residual_norm_max_uT"],
                      "field_norm":[mag["field_norm_min_uT"],mag["field_norm_max_uT"]]},
  "measurement_model":{"accel_std_max":C["sensor_model"]["accel_effective_std_max_mps2"]}}
