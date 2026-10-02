"""Norm envelope for the non-INJ covariance remainder.

Write P=P_struct+E, ||E||_2<=D.  This module propagates D without pretending
that reset-generated perpendicular directions remain in the three-scalar INJ
algebra.  Bounds are subordinate to the exact structured/vector path.
"""
from __future__ import annotations
import math

def prediction_defect(D,F_norm,axis_rebase_defect=0.,process_defect=0.):
 """P -> FPF'+Q, plus the separately derived changing-axis rebase remainder."""
 return F_norm*F_norm*D+axis_rebase_defect+process_defect

def correction_defect(D,P_struct_norm,H_norm,S_struct_inv_norm):
 """Resolvent bound for C(P)=P-PH'(HPH'+R)^-1HP.

 The structured innovation inverse is S0^-1.  If
 q=||S0^-1|| ||H||^2 D < 1, Banach's lemma bounds the perturbed inverse.
 This is a sufficient perturbation envelope, not a claim of losslessness.
 """
 if D<0:raise ValueError("D must be nonnegative")
 q=S_struct_inv_norm*H_norm*H_norm*D
 if q>=1:return {"closed":False,"D_next":math.inf,"q":q,
                 "reason":"innovation resolvent radius exhausted"}
 si2=S_struct_inv_norm/(1-q)
 invdiff=S_struct_inv_norm*H_norm*H_norm*D*si2
 p=P_struct_norm;e=D;h=H_norm
 # Difference of P H' S^-1 H P, retaining the three product differences.
 dt=(e*h*h*si2*(p+e)
     +p*h*h*invdiff*(p+e)
     +p*h*h*S_struct_inv_norm*e)
 return {"closed":True,"D_next":D+dt,"q":q,
         "perturbed_inverse_norm":si2}

def reset_defect(D,G_norm,new_structured_reset_defect):
 """Congruence reset plus exact structured-reset split remainder."""
 return G_norm*G_norm*D+new_structured_reset_defect

def chronology_step(D,operations):
 """Apply precomputed literal-order defect operations.

 operation dictionaries:
 prediction: F_norm, axis_rebase_defect, process_defect
 correction: P_struct_norm,H_norm,S_struct_inv_norm
 reset: G_norm,new_structured_reset_defect
 """
 z=float(D);trace=[]
 for op in operations:
  kind=op["kind"]
  if kind=="prediction":
   z=prediction_defect(z,op["F_norm"],op.get("axis_rebase_defect",0.),
                       op.get("process_defect",0.));row={"closed":True,"D_next":z}
  elif kind=="correction":
   row=correction_defect(z,op["P_struct_norm"],op["H_norm"],op["S_struct_inv_norm"])
   z=row["D_next"]
  elif kind=="reset":
   z=reset_defect(z,op["G_norm"],op["new_structured_reset_defect"])
   row={"closed":True,"D_next":z}
  else:raise ValueError("unknown covariance operation "+str(kind))
  trace.append({"kind":kind,**row})
  if not row["closed"]:break
 return {"D_next":z,"trace":trace,"closed":all(x["closed"] for x in trace),
         "classification_if_exhausted":"D_SUFFICIENT_BOUND_FAILURE"}
