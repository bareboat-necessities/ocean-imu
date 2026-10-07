"""Conditional parity-block information ceilings, with explicit relaxations.

Rows are assembled in the exact 21-state coordinates then restricted to the
exact E/O parity blocks. The bounds dominate H' R_floor^-1 H. Larger R only decreases information.
The accelerometer bound is a block-diagonal Cauchy relaxation, NOT an actual
measurement row and NOT a claim that independent block extrema are reachable.

Planar assumptions used here: lever arm disabled in the shipping/default
history, BA active after refinement, |f_cog|<=30 candidate bound, Racc std
>=0.2 on the replay/default profile, Rmag std=0.8, vertical R_S std>=0.15 with literal axis factors (0.72,0.50,1).
The force and diagnostic noise-profile bindings remain unproved for all time.
"""
from __future__ import annotations
import numpy as np
from .planar_parity import EVEN,ODD

def skew(v):
 x,y,z=v;return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def restrict(J,idx): return J[np.ix_(idx,idx)]

def acc_ceiling(*,f_norm=30.,racc_std=.2,use_ba=True):
 # For A_i with ||A_i||<=w_i, ||sum A_i x_i||^2
 # <= (sum w_i) sum w_i ||x_i||^2. An outer product w w' tensor I
 # is NOT a bound when the A_i have different orientations/nullspaces.
 if not np.isfinite(f_norm) or f_norm<0 or not np.isfinite(racc_std) or racc_std<=0:
  raise ValueError("finite nonnegative force and positive noise required")
 blocks=[(0,3,f_norm),(15,18,1.)]
 if use_ba: blocks.append((18,21,1.))
 J=np.zeros((21,21));total=sum(w for _,_,w in blocks)
 for a,b,w in blocks:J[a:b,a:b]=total*w/(racc_std**2)*np.eye(3)
 return restrict(J,EVEN),restrict(J,ODD)

def S_ceiling(*,rs_std_min=.15,axis_factors=(.72,.50,1.)):
 std=rs_std_min*np.asarray(axis_factors,dtype=float)
 if std.shape!=(3,) or not np.isfinite(std).all() or np.any(std<=0):
  raise ValueError("three positive finite S noise floors required")
 J=np.zeros((21,21));J[12:15,12:15]=np.diag(1/std**2)
 return restrict(J,EVEN),restrict(J,ODD)

def mag_ceiling(*,B=75.,rmag_std=.8):
 # H=-[B]x, H'H <= B^2 I in attitude block.
 J=np.zeros((21,21));J[:3,:3]=B*B/(rmag_std**2)*np.eye(3)
 return restrict(J,EVEN),restrict(J,ODD)

def certificate():
 ae,ao=acc_ceiling();se,so=S_ceiling();me,mo=mag_ceiling()
 return {"qualification":"OU3_PLANAR_INFORMATION_CEILINGS_V2",
         "assumptions":{"f_cog_norm_upper":30.0,"Racc_std_lower":.2,"Rmag_std":.8,"R_S_vertical_std_lower":.15,"R_S_axis_factors":[.72,.50,1.],
                        "lever_arm_disabled":True,"BA_active":True},
         "even":{"acc_lambda_max":float(np.linalg.eigvalsh(ae).max()),"S_lambda_max":float(np.linalg.eigvalsh(se).max()),"mag_lambda_max":float(np.linalg.eigvalsh(me).max())},
         "odd":{"acc_lambda_max":float(np.linalg.eigvalsh(ao).max()),"S_lambda_max":float(np.linalg.eigvalsh(so).max()),"mag_lambda_max":float(np.linalg.eigvalsh(mo).max())},
         "matrix_ceiling_method":"weighted block-diagonal Cauchy majorant; cross correlation relaxed explicitly",
         "all_time_f_bound_verified":False,
         "literal_profile_binding_verified":False,
         "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
