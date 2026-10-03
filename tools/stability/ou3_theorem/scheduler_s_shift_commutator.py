"""Local covariance sensitivity to shifting one S=0 correction by one sample.

For maps prediction A(P)=FPF'+Q and S correction R(P), neighboring scheduler
words replace ...R∘A... by ...A∘R... locally. This module computes the exact
commutator on a covariance and provides spectral norm used for branch-tube
radius propagation. Rigorous promotion requires intervalizing over the local
phase covariance tube.
"""
from __future__ import annotations
import numpy as np
from .planar_periodic_covariance_tube import predict,joseph
def shift_defect(P,F,Q,H,R):
 a=joseph(predict(P,F,Q),H,R)
 b=predict(joseph(P,H,R),F,Q)
 return a-b
def defect_norm(P,F,Q,H,R):
 D=shift_defect(P,F,Q,H,R);return float(np.linalg.norm((D+D.T)/2,2))
def certificate():
 return {"qualification":"OU3_SCHEDULER_S_SHIFT_COMMUTATOR_V1",
         "identity":"neighboring scheduler branches differ locally by exchanging one S correction with one prediction/correction sample boundary",
         "adjacent_words_hamming_at_most_two":True,
         "interval_commutator_bound_verified":False,"theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
