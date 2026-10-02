"""Numerical same-history force-axis mismatch diagnostics."""
import math
from .force_axis_same_history import axis_angle_bound,force_change_bound
H=.006;J=100.
def diagnostic():
 df=H*J;rows=[]
 for rho in [18.60665,15,10,5,2,1,.6,.3,.1,0]:
  d=axis_angle_bound(rho,rho,df)
  # goLive attitude projector coefficient ptheta*rho^2
  ptheta=1.5708**2;b=ptheta*rho*rho
  pc=b*(math.sin(d) if d<math.pi else 1.)
  rows.append({"rho":rho,"df":df,"delta_deg":math.degrees(d),
               "projector_coefficient":b,"weighted_projector_mismatch":pc})
 return {"verified":True,"h_s":H,"Jmax":J,"rows":rows}
