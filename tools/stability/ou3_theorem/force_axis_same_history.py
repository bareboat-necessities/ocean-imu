"""Same-history force-axis mismatch bound with nonsingular low-force treatment."""
from __future__ import annotations
import math

def axis_angle_bound(rho0,rho1,df):
 """Sharp chord/angle bound from two nonzero vectors with ||f1-f0||<=df.

 If r=min(rho0,rho1)>0 and df<2r, angle <=2 asin(df/(2r)).
 """
 r=min(rho0,rho1)
 if r<=0:return math.pi
 if df>=2*r:return math.pi
 return 2*math.asin(min(1.,df/(2*r)))

def force_change_bound(h,Jmax,frame_rate_term=0.):
 """Same-history force-vector change.

 World physical acceleration contributes h Jmax.  If comparing in a rotating
 local frame, a separately derived frame-rate contribution may be added as
 h*Omega*rho; callers must not independently box it.
 """
 return h*Jmax+frame_rate_term

def weighted_rebase_bound(b,c,rho0,rho1,h,Jmax):
 """No singularity at rho->0: anisotropic projector coefficient b scales rho^2.

 For A=aI+b nn'+c[n]x, use exact angle when resolvable; otherwise worst basis
 difference. This remains finite as force vanishes.
 """
 df=h*Jmax;delta=axis_angle_bound(rho0,rho1,df)
 return {"delta_rad":delta,"delta_deg":math.degrees(delta),
  "projector_cost":abs(b)*math.sin(delta) if delta<math.pi else abs(b),
  "skew_cost":2*abs(c)*math.sin(delta/2),
  "low_force_safe":True,"df":df}

def threshold_split(h,Jmax,eps):
 """rho threshold above which angular mismatch <=eps from jerk alone."""
 s=math.sin(eps/2)
 return math.inf if s<=0 else h*Jmax/(2*s)
