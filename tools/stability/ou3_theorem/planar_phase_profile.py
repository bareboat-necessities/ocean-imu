"""Generate compact phase profile from a literal replay export.

A rigorous tube certificate needs the actual phase-indexed covariance centers.
This module validates a profile JSON emitted by the native probe and computes
candidate Loewner radii from adjacent source periods. It never promotes replay
radii to a theorem.
"""
from __future__ import annotations
import json,numpy as np
def compare_profiles(a,b):
    worst=0.
    for A,B in zip(a,b):
        D=np.asarray(A)-np.asarray(B);worst=max(worst,float(np.linalg.norm(D,2)))
    return worst
def certificate():
 return {"qualification":"OU3_PLANAR_PHASE_PROFILE_V1","profile_available":False,
         "candidate_radius_method":"spectral difference between consecutive settled 20-s source profiles, inflated before interval containment",
         "replay_radius_is_not_certificate":True,"theorem_closed":False}
if __name__=="__main__":print(json.dumps(certificate(),indent=2,sort_keys=True))
