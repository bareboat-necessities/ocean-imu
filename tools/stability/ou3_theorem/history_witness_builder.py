"""Build nonlinear and kernel witnesses from one propagated history cell.

No midpoint witness admission.  All outputs retain the root dependency token.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .marine_magnetic_qcqp import VectorBox,MagneticEventBox,MarineMagneticWitness
from .compatibility_line_interval import joint_line,line_angle_radius

@dataclass(frozen=True)
class TimedVectorBox:
 t_lo:float;t_hi:float;value:VectorBox

@dataclass(frozen=True)
class AttitudeGravityBox:
 t_lo:float;t_hi:float;gravity:VectorBox

@dataclass(frozen=True)
class KernelLineCone:
 midpoint:np.ndarray;component_radius:np.ndarray;angle_rad:float;verified:bool

@dataclass
class HistoryWitness:
 dependency_token:str
 marine_magnetic:MarineMagneticWitness
 kernel_line:KernelLineCone

def _covering_pairs(samples,T):
 """All endpoint pairs whose time separation certainly spans a complete T window.

 A real MARINE diameter witness may occur at interior times. The causal
 propagator must therefore tag candidate witness samples densely enough; this
 routine returns every certainly separated pair rather than choosing a midpoint.
 """
 out=[]
 for i,a in enumerate(samples):
  for b in samples[i+1:]:
   if b.t_lo-a.t_hi>=T-1e-12:out.append((a.value,b.value))
 return out

def build_history_witness(*,dependency_token,velocity,position,acceleration,jerk,
                          gravity,magnetic_events,start,end,
                          attitude_line_mid,attitude_line_rad,ba_graph_mid,ba_graph_rad,
                          marine_T=30.):
 if not dependency_token:raise ValueError("dependency token required")
 seqs=(velocity,position,acceleration,jerk,gravity)
 if any(not x for x in seqs):raise ArithmeticError("causal physical witness history incomplete")
 gp=_covering_pairs(gravity,marine_T);pp=_covering_pairs(position,marine_T)
 if end-start>=marine_T and (not gp or not pp):raise ArithmeticError("30-s witness candidates missing")
 # Joint graph line interval: r=(a,-G a). Midpoint and first-order interval
 # radius with rigorous product remainder |dG| |da|.
 a=np.asarray(attitude_line_mid,float);ar=np.asarray(attitude_line_rad,float)
 G=np.asarray(ba_graph_mid,float);Gr=np.asarray(ba_graph_rad,float)
 r=joint_line(a,G); rr=np.zeros(21);rr[:3]=ar
 rr[18:21]=np.abs(G)@ar+Gr@np.abs(a)+Gr@ar
 ac=line_angle_radius(r,rr)
 if not ac["verified"]:raise ArithmeticError("compatibility line cone includes zero")
 mm=MarineMagneticWitness([x.value for x in velocity],[x.value for x in position],
      [x.value for x in acceleration],[x.value for x in jerk],gp,pp,list(magnetic_events),start,end)
 return HistoryWitness(dependency_token,mm,KernelLineCone(r,rr,ac["angle_rad"],True))
