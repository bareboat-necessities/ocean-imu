"""Literal parity-block correction information ceilings for planar history.

Rows are assembled in the exact 21-state coordinates then restricted to the
exact E/O parity blocks. The ceiling is Jbar=H' R_floor^-1 H. It is an upper
information matrix: larger R only decreases information.

Planar assumptions used here: lever arm disabled in the shipping/default
history, BA active after refinement, |f_cog|<=30 candidate bound, Racc std
>=0.2 on the replay/default profile, Rmag std=0.8, R_S std>=0.15 with axis
factors already included by the committed R_S lower clamp.
"""
from __future__ import annotations
import numpy as np
from .planar_parity import EVEN,ODD

def skew(v):
 x,y,z=v;return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def restrict(J,idx): return J[np.ix_(idx,idx)]

def acc_ceiling(*,f_norm=30.,racc_std=.2,use_ba=True):
 # Worst orientation-free Loewner ceiling for Jatt=-[f]x, Jaw=R, Jba=I.
 # Build block row norm by replacing actual Jatt with f_norm*orthogonal action.
 # H'H <= (||Jatt||+||Jaw||+||Jba||)^2 on the coupled row; retain state blocks
 # via exact Cauchy outer bound vv' tensor I, which keeps cross terms.
 blocks=[(0,3,f_norm),(15,18,1.)]
 if use_ba: blocks.append((18,21,1.))
 J=np.zeros((21,21)); r=racc_std**2
 for a,b,wa in blocks:
  for c,d,wb in blocks:J[a:b,c:d]+=wa*wb/r*np.eye(3)
 return restrict(J,EVEN),restrict(J,ODD)

def S_ceiling(*,rs_std_min=.15):
 J=np.zeros((21,21));J[12:15,12:15]=np.eye(3)/(rs_std_min**2)
 return restrict(J,EVEN),restrict(J,ODD)

def mag_ceiling(*,B=75.,rmag_std=.8):
 # H=-[B]x, H'H <= B^2 I in attitude block.
 J=np.zeros((21,21));J[:3,:3]=B*B/(rmag_std**2)*np.eye(3)
 return restrict(J,EVEN),restrict(J,ODD)

def certificate():
 ae,ao=acc_ceiling();se,so=S_ceiling();me,mo=mag_ceiling()
 return {"qualification":"OU3_PLANAR_LITERAL_INFORMATION_CEILINGS_V1",
         "assumptions":{"f_cog_norm_upper":30.0,"Racc_std_lower":.2,"Rmag_std":.8,"R_S_std_lower":.15,
                        "lever_arm_disabled":True,"BA_active":True},
         "even":{"acc_lambda_max":float(np.linalg.eigvalsh(ae).max()),"S_lambda_max":float(np.linalg.eigvalsh(se).max()),"mag_lambda_max":float(np.linalg.eigvalsh(me).max())},
         "odd":{"acc_lambda_max":float(np.linalg.eigvalsh(ao).max()),"S_lambda_max":float(np.linalg.eigvalsh(so).max()),"mag_lambda_max":float(np.linalg.eigvalsh(mo).max())},
         "matrix_ceiling_method":"coupled block Cauchy outer product; cross terms retained",
         "all_time_f_bound_verified":False,
         "literal_profile_binding_verified":False,
         "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
