"""Outer interpolation of shared history knots to delivered-sample cells."""
from __future__ import annotations
import re,math
from .causal_tuner_interval import I

_PAT=re.compile(r"^(.*)_([0-2])_([0-9]+(?:\.[0-9]+)?)$")

def parse_knots(cell):
 out={}
 for name,x in cell.coordinates:
  m=_PAT.match(name)
  if not m:continue
  kind,axis,t=m.group(1),int(m.group(2)),float(m.group(3))
  out.setdefault((kind,axis),[]).append((t,I(x.lo,x.hi)))
 for k in out:out[k].sort()
 return out

def lipschitz_outer(knots,t,L,amplitude=None):
 """Outer value at t from ALL knot constraints |x(t)-x(ti)|<=L|dt|."""
 candidates=[]
 for ti,x in knots:
  candidates.append(I(x.lo-L*abs(t-ti),x.hi+L*abs(t-ti)))
 lo=max(x.lo for x in candidates);hi=min(x.hi for x in candidates)
 if amplitude is not None:lo=max(lo,-amplitude);hi=min(hi,amplitude)
 if lo>hi:raise ArithmeticError("knot cell violates Lipschitz history")
 return I(lo,hi)

def bounded_outer(knots,t,amplitude):
 # Without derivative information, knot interpolation cannot reduce the theorem
 # amplitude between knots. Returning the full envelope is rigorous.
 return I(-amplitude,amplitude)

def sample_history(cell,times,constants):
 k=parse_knots(cell);imu=constants["imu_bias"];m=constants["marine_motion"];out=[]
 for t in times:
  row={"t":t}
  for a in range(3):
   row[f"slow_accel_{a}"]=lipschitz_outer(k[("slow_accel",a)],t,imu["D_a_s_mps3"],imu["B_a_s_mps2"])
   row[f"slow_gyro_{a}"]=lipschitz_outer(k[("slow_gyro",a)],t,imu["D_g_s_rad_s2"],imu["B_g_s_rad_s"])
   # Physical p/v are cross-constrained by derivatives. Use both direct knot
   # Lipschitz envelopes; later causal consistency intersects p'=v, v'=a.
   row[f"physical_p_{a}"]=lipschitz_outer(k[("physical_p",a)],t,m["V_max_mps"],m["P_max_m"])
   row[f"physical_v_{a}"]=lipschitz_outer(k[("physical_v",a)],t,m["A_max_mps2"],m["V_max_mps"])
   row[f"gravity_dir_{a}"]=bounded_outer(k[("gravity_dir",a)],t,1.)
   row[f"mag_field_{a}"]=bounded_outer(k[("mag_field",a)],t,constants["magnetic_service"]["field_norm_max_uT"])
  out.append(row)
 return out

def fast_primitive_increment_outer(cell,kind,axis,t0,t1,constants):
 """Use literal all-window cap, not interpolation of primitive knots."""
 imu=constants["imu_bias"];T=t1-t0
 if kind=="accel":B=imu["B_a_f_mps2"];H=imu["fast_accel_horizon_s"];C=imu["fast_accel_accumulation_cap_mps"]
 else:B=imu["B_g_f_rad_s"];H=imu["fast_gyro_horizon_s"];C=imu["fast_gyro_accumulation_cap_rad"]
 if T<0 or T>H+1e-12:raise ArithmeticError("FAST increment outside qualified horizon")
 cap=min(B*T,C);return I(-cap,cap)

def acceleration_outer_from_velocity_knots(cell,axis,t,constants,probe_s=None):
 """Causal outer acceleration from the SAME physical-v history and jerk bound.

 For any h>0, v(t+h)-v(t)=integral a and |a(s)-a(t)|<=J|s-t| imply
 |a(t)-(v(t+h)-v(t))/h|<=J h/2.  Intersect forward/backward knot pairs.
 This retains signed primitive information that a pointwise A_max cube loses.
 """
 k=parse_knots(cell);vk=k[("physical_v",axis)]
 A=constants["marine_motion"]["A_max_mps2"];J=constants["marine_motion"]["J_max_mps3"]
 cand=[I(-A,A)]
 for ti,vi in vk:
  h=ti-t
  if abs(h)<1e-12:continue
  # Bound v(t) from all velocity knots before forming the signed secant.
  vt=lipschitz_outer(vk,t,A,constants["marine_motion"]["V_max_mps"])
  if h>0:
   sec=I((vi.lo-vt.hi)/h,(vi.hi-vt.lo)/h)
  else:
   hh=-h;sec=I((vt.lo-vi.hi)/hh,(vt.hi-vi.lo)/hh)
  rad=.5*J*abs(h);cand.append(I(sec.lo-rad,sec.hi+rad))
 lo=max(x.lo for x in cand);hi=min(x.hi for x in cand)
 if lo>hi:raise ArithmeticError("velocity knots incompatible with jerk/acceleration history")
 return I(lo,hi)

def physical_acceleration_history(cell,times,constants):
 return [tuple(acceleration_outer_from_velocity_knots(cell,a,t,constants) for a in range(3))
         for t in times]
