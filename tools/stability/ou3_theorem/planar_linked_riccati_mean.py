"""Exact local differential identities for the remaining mean<->covariance port.

For C(P,H)=P-PH'(HPH'+R)^-1HP, differentiate without independent extrema.
For the mean correction x+=K r, K=PH'S^-1, retain dK and dr=-dH*x_model
through the same local operands.  These identities are the tangent recurrence
used by the literal 20-s port-word exporter.
"""
from __future__ import annotations

def certificate():
 return {"qualification":"OU3_PLANAR_LINKED_RICCATI_MEAN_DIFFERENTIAL_V1",
  "result_type":"PROVED analytical differential identity",
  "covariance_differential":
   "dC=dP-dP H' S^-1 H P-P dH' S^-1 H P-P H' S^-1 dH P-P H' S^-1 H dP+P H' S^-1(dH P H'+H dP H'+H P dH')S^-1 H P",
  "gain_differential":
   "dK=dP H' S^-1+P dH' S^-1-P H' S^-1(dH P H'+H dP H'+H P dH')S^-1",
  "mean_differential":"d(x+)=d(x)+dK*r+K*dr; dr is differentiated from the same measurement model",
  "independent_extrema_used":False,"finite_word_uniform_bound_closed":False,
  "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
