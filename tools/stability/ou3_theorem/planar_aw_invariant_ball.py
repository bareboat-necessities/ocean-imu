"""Planar AW invariant-ball sufficient condition.

For accepted accel updates aw+ = phi aw + Kaw r, with
r=f_meas-R(aw-g)-ba on the planar no-lever history. A norm-only sufficient
invariant radius A is
 phi*A + k*(F + A + g + B) <= A,
i.e. A >= k*(F+g+B)/(1-phi-k), requiring phi+k<1.
This is deliberately conservative; if it fails, use the signed matrix
I-Kaw R rather than treating feedback as forcing.
"""
from __future__ import annotations
def sufficient_radius(phi,k,F,g=9.80665,B=.2):
 if phi+k>=1:return None
 return k*(F+g+B)/(1-phi-k)
def certificate():
 return {"qualification":"OU3_PLANAR_AW_INVARIANT_BALL_V1",
         "norm_only_condition":"phi+k<1",
         "warning":"failure of norm-only condition is D-type; exact update contains stabilizing -Kaw R aw feedback",
         "signed_matrix_certificate_verified":False,"theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
