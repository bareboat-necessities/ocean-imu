# ruff: noqa: F401, F811
"""Causal history-cell witness builder joining nonlinear premises and shaped ratio."""
from __future__ import annotations
from dataclasses import dataclass,field
import numpy as np
from .reachable_history_enclosure import HistoryCell
from .marine_magnetic_qcqp import VectorBox,MagneticEventBox,MarineMagneticWitness,callback
from .compatibility_line_interval import joint_line,line_angle_radius
from .constructive_cell_ratio import kernel_certificate_from_line_cone,combine
from .verified_linked_qcqp import Box,verified_interval_supply_bnb
from .temporal_source_domain import build_temporal_domain,nonlinear_side_constraints

@dataclass
class CausalWitnessState:
 root:HistoryCell
 times:list[float]=field(default_factory=list)
 velocity:list[VectorBox]=field(default_factory=list)
 position:list[VectorBox]=field(default_factory=list)
 acceleration:list[VectorBox]=field(default_factory=list)
 jerk:list[VectorBox]=field(default_factory=list)
 gravity:list[VectorBox]=field(default_factory=list)
 magnetic_events:list[MagneticEventBox]=field(default_factory=list)
 attitude_lines:list[np.ndarray]=field(default_factory=list)
 ba_graphs:list[np.ndarray]=field(default_factory=list)
 line_radii:list[np.ndarray]=field(default_factory=list)

 def append_physical(self,t,v,p,a,j,g):
  if self.times and t<=self.times[-1]:raise ValueError("strict causal time required")
  self.times.append(float(t));self.velocity.append(v);self.position.append(p);self.acceleration.append(a);self.jerk.append(j);self.gravity.append(g)

 def append_magnetic(self,event):
  if event.time_lo<0 or event.time_hi<event.time_lo:raise ValueError("magnetic event time")
  self.magnetic_events.append(event)

 def append_kernel_line(self,attitude_line,ba_graph,component_radius):
  self.attitude_lines.append(np.asarray(attitude_line,float));self.ba_graphs.append(np.asarray(ba_graph,float));self.line_radii.append(np.asarray(component_radius,float))

 def _window_pairs(self,T,values,kind,minimum):
  pairs=[]
  if not self.times:return pairs
  # For each complete T window rooted at a sampled time, require a certified
  # pair inside that window. Search all endpoint pairs and retain one only when
  # its whole-box lower span reaches the premise. Absence is represented by a
  # degenerate unresolved pair so callback cannot silently pass.
  from .marine_magnetic_qcqp import difference_span_witness,gravity_direction_span_witness,TRUE
  for i,t0 in enumerate(self.times):
   if t0+T>self.times[-1]+1e-12:break
   ids=[j for j,t in enumerate(self.times) if t0-1e-12<=t<=t0+T+1e-12]
   found=None
   for a in ids:
    for b in ids:
     q=gravity_direction_span_witness(values[a],values[b],minimum) if kind=="gravity" else difference_span_witness(values[a],values[b],minimum)
     if q is TRUE:found=(values[a],values[b]);break
    if found:break
   if found is None:
    # Callback returns UNKNOWN/FALSE for identical broad boxes; never fabricate excitation.
    found=(values[ids[0]],values[ids[0]])
   pairs.append(found)
  return pairs

 def marine_magnetic_witness(self,constants):
  if len(self.times)<2:raise ValueError("physical history absent")
  m=constants["marine_motion"]
  gp=self._window_pairs(m["attitude_excitation"]["T_E_s"],self.gravity,"gravity",m["attitude_excitation"]["theta_E_rad"])
  pp=self._window_pairs(m["displacement_excitation"]["T_P_s"],self.position,"position",m["displacement_excitation"]["P_E_m"])
  return MarineMagneticWitness(self.velocity,self.position,self.acceleration,self.jerk,gp,pp,self.magnetic_events,self.times[0],self.times[-1])

 def compatibility_cone(self,index=-1):
  r=joint_line(self.attitude_lines[index],self.ba_graphs[index])
  rad=np.zeros(21);rad[:3]=self.line_radii[index][:3]
  # BA line uncertainty from A_ba*a is supplied by builder in final three entries if present.
  if len(self.line_radii[index])>=6:rad[18:21]=self.line_radii[index][3:6]
  c=line_angle_radius(r,rad)
  if not c["verified"]:raise ArithmeticError("compatibility cone not separated")
  return r,rad,c

def certify_leaf(state,constants,iq,A0,T,A1,line1,e_box,u_box,symbols,source_side_builder,max_leaves=200000):
 """Join BOTH certificate halves on one dependency token."""
 premise=callback(state.marine_magnetic_witness(constants),constants)
 if premise is not True:
  return {"verified":False,"reason":"nonlinear premise unresolved" if premise is None else "premise false"}
 line0,rad,cone=state.compatibility_cone()
 kc=kernel_certificate_from_line_cone(A0,line0,rad,T,A1,line1)
 if not kc["zero_set_excluded"]:return {"verified":False,"reason":"kernel restricted action"}
 domain=build_temporal_domain(symbols,state.times[-1]-state.times[0]);LA,Lb=domain.matrices()
 side=source_side_builder(state,constants)
 supply=verified_interval_supply_bnb(iq.B,iq.C,e_box,u_box,LA,Lb,side,max_leaves=max_leaves)
 ratio=combine(state.root.prefix_token,iq,kc,supply)
 return {"verified":ratio.verified,"ratio":ratio,"kernel":kc,"supply":supply,
         "same_dependency_token":iq.dependency_token==state.root.prefix_token}
