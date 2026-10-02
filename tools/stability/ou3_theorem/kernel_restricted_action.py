"""Interval restriction by a word-dependent compatibility line.

The line is supplied by the literal complete-word kernel construction for the
SAME history cell.  No fixed global vector is assumed.
"""
from __future__ import annotations
import math
import numpy as np
from .interval_riccati_21 import IMat,matmul,transpose
from .rank_loss_interval_factor import exact
from .interval_riccati import symmetric_interval_gershgorin

def householder_complement(r):
 r=np.asarray(r,float).reshape(-1);n=np.linalg.norm(r)
 if not(n>0 and np.isfinite(n)):raise ValueError("nonzero finite kernel line")
 r=r/n
 # Deterministic Householder sends r to +/- e_j; remaining columns span r^perp.
 j=int(np.argmax(np.abs(r)));sgn=1. if r[j]>=0 else -1.
 v=r.copy();v[j]+=sgn;nv=np.linalg.norm(v)
 if nv<1e-14:
  Q=np.eye(len(r))
 else:
  v/=nv;Q=np.eye(len(r))-2*np.outer(v,v)
 keep=[k for k in range(len(r)) if k!=j]
 U=Q[:,keep]
 if np.linalg.norm(U.T@r)>2e-12:raise ArithmeticError("complement construction")
 return U

def restrict_interval_action(A:IMat,line_mid,line_radius=None):
 """Restrict A to a certified complement of an interval kernel line.

 Nonzero line radius requires a separate angle certificate; fail closed rather
 than pretending the midpoint complement is valid for the whole cell.
 """
 r=np.asarray(line_mid,float)
 if line_radius is not None and np.max(np.abs(line_radius))>0:
  from .compatibility_line_interval import line_angle_radius
  from .angle_projector_enclosure import projector_cone_restricted_lower
  ac=line_angle_radius(r,line_radius)
  if not ac["verified"]: raise ArithmeticError("kernel line cone not separated from zero")
  c=projector_cone_restricted_lower(A,r,ac["angle_rad"])
  if not c["verified"]: raise ArithmeticError("angle-aware restricted action not positive")
  U=householder_complement(r);Ui=exact(U.tolist())
  R=matmul(matmul(transpose(Ui),A),Ui)
  return R,{**c,"kernel_dimension":1,"line_component_radius_used":True}
 U=householder_complement(r);Ui=exact(U.tolist())
 R=matmul(matmul(transpose(Ui),A),Ui)
 lo,hi=symmetric_interval_gershgorin(R.mid,R.rad)
 return R,{"verified":lo>0,"restricted_eigen_lower":max(0.,lo),
           "restricted_eigen_upper":hi,"kernel_dimension":1,
           "fixed_global_kernel_used":False}

def augment_later_word_action(A0:IMat,T:IMat,A1:IMat):
 """A_super=A0+T' A1 T in root coordinates."""
 from .interval_riccati_21 import add
 return add(A0,matmul(matmul(transpose(T),A1),T))

def certify_rotating_line(A0:IMat,line0,T:IMat,A1:IMat,line1):
 """Use later-word action before restricting the first compatibility line."""
 As=augment_later_word_action(A0,T,A1)
 R,c=restrict_interval_action(As,line0)
 # Also verify the carried first-line image is not silently declared equal to
 # the next line. Quantitative separation is reported, not assumed.
 tm=np.asarray(T.mid,float)@np.asarray(line0,float);r1=np.asarray(line1,float)
 nt=np.linalg.norm(tm);nr=np.linalg.norm(r1)
 if nt>0 and nr>0:
  cos=abs(float(tm@r1)/(nt*nr));cos=min(1.,max(0.,cos));ang=math.acos(cos)
 else:ang=0.
 c["midpoint_line_rotation_rad"]=ang;c["later_word_action_included"]=True
 return R,c
