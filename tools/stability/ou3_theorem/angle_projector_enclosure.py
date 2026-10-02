"""Angle-aware restriction of an interval action to a line-complement cone."""
from __future__ import annotations
import math
import numpy as np
from .interval_riccati_21 import IMat
from .kernel_restricted_action import householder_complement
from .rank_loss_interval_factor import exact
from .interval_riccati_21 import matmul,transpose
from .interval_riccati import symmetric_interval_gershgorin

def interval_spectral_norm_upper(A:IMat):
 # ||A||2 <= ||mid||2 + ||rad||2 <= ||mid||2 + sqrt(||rad||1||rad||inf)
 M=np.asarray(A.mid,float);R=np.asarray(A.rad,float)
 return float(np.linalg.norm(M,2)+math.sqrt(np.linalg.norm(R,1)*np.linalg.norm(R,np.inf)))

def projector_cone_restricted_lower(A:IMat,line_mid,angle_rad:float):
 """Uniform lower bound on x'Ax for unit x perpendicular to any line in cone.

 Let P0=I-r0r0'. For unit lines r with angle<=delta,
 ||P-P0||2=sin(delta). For x in range(P), y=P0 x has
 ||x-y||<=eps and ||y||>=sqrt(1-eps^2). Compare quadratic forms against
 the midpoint-complement restricted floor and full action norm.
 """
 if not(math.isfinite(angle_rad) and 0<=angle_rad<math.pi/2):return {"verified":False,"reason":"invalid line cone"}
 r=np.asarray(line_mid,float);r/=np.linalg.norm(r)
 U=householder_complement(r);Ui=exact(U.tolist());R0=matmul(matmul(transpose(Ui),A),Ui)
 lam0,_=symmetric_interval_gershgorin(R0.mid,R0.rad)
 normA=interval_spectral_norm_upper(A);eps=math.sin(angle_rad)
 # x=y+z, ||z||<=eps, ||y||>=sqrt(1-eps^2).
 # x'Ax >= lam0||y||^2 - 2||A||||y||||z||-||A||||z||^2.
 # Using ||y||<=1 in the adverse cross term gives rigorous simple floor.
 lower=lam0*(1-eps*eps)-normA*(2*eps+eps*eps)
 return {"verified":lower>0,"restricted_eigen_lower":max(0.,math.nextafter(lower,-math.inf)),
         "midpoint_restricted_lower":lam0,"action_norm_upper":normA,
         "line_angle_rad":angle_rad,"projector_distance_upper":eps,
         "fixed_global_kernel_used":False,"whole_line_cone_covered":True}
