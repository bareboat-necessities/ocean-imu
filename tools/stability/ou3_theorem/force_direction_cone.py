"""Magnitude + directional-cone force geometry for accelerometer proof."""
from __future__ import annotations
from dataclasses import dataclass
import math,numpy as np
from .interval_riccati_21 import IMat,N

def skew(x):
 x=np.asarray(x,float);return np.array([[0,-x[2],x[1]],[x[2],0,-x[0]],[-x[1],x[0],0.]])

@dataclass(frozen=True)
class ForceCone:
 rho_lo:float;rho_hi:float
 axis:np.ndarray
 half_angle:float
 def __post_init__(self):
  if not(0<=self.rho_lo<=self.rho_hi):raise ValueError("force magnitude interval")
  if not(0<=self.half_angle<math.pi/2):raise ValueError("direction cone")
  if abs(np.linalg.norm(self.axis)-1)>1e-10:raise ValueError("unit cone axis")
 def midpoint_force(self):
  return .5*(self.rho_lo+self.rho_hi)*self.axis
 def skew_operator_radius(self):
  """Exact operator perturbation bound about midpoint rho*n0.

  |rho-rhoc| contribution <= drho.
  For unit n within angle alpha, ||n-n0||<=2 sin(alpha/2), hence
  ||[rho n-rhoc n0]x||2 = ||rho n-rhoc n0||
     <= drho + rho_hi*2 sin(alpha/2).
  """
  dr=.5*(self.rho_hi-self.rho_lo)
  return dr+self.rho_hi*2*math.sin(.5*self.half_angle)

def accel_h_force_cone(cone:ForceCone,R_mid,R_rad):
 hm=np.zeros((3,N));hr=np.zeros((3,N))
 f0=cone.midpoint_force();hm[:,:3]=-skew(f0)
 # Convert spectral skew perturbation to safe entrywise radius.
 hr[:,:3]=cone.skew_operator_radius()
 hm[:,15:18]=np.asarray(R_mid,float);hr[:,15:18]=np.asarray(R_rad,float)
 for i in range(3):hm[i,18+i]=1.
 return IMat(tuple(map(tuple,hm)),tuple(map(tuple,hr)))

_OCTAHEDRAL_HALF_ANGLE=math.acos(1/math.sqrt(3))
def octahedral_direction_cones(half_angle=_OCTAHEDRAL_HALF_ANGLE):
 """Six axis cones cover S^2 when half-angle >= acos(1/sqrt(3))."""
 axes=[]
 for i in range(3):
  for s in (-1.,1.):
   a=np.zeros(3);a[i]=s;axes.append(a)
 return [ForceCone(0.,18.60665,a,half_angle) for a in axes]

def split_radial(c:ForceCone):
 m=.5*(c.rho_lo+c.rho_hi)
 return ForceCone(c.rho_lo,m,c.axis,c.half_angle),ForceCone(m,c.rho_hi,c.axis,c.half_angle)
