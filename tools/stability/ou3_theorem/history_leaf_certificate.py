"""Connect one causal history witness to nonlinear and kernel certificates."""
from __future__ import annotations
from .marine_magnetic_qcqp import callback
from .kernel_restricted_action import augment_later_word_action,restrict_interval_action
from .constructive_cell_ratio import combine

def certify_leaf(*,witness,constants,A0,T,A1,interval_joint,supply_certificate):
 side=callback(witness.marine_magnetic,constants)
 if side is not True:
  return {"verified":False,"reason":"nonlinear premises unresolved" if side is None else "outside nonlinear premises"}
 As=augment_later_word_action(A0,T,A1)
 R,k=restrict_interval_action(As,witness.kernel_line.midpoint,witness.kernel_line.component_radius)
 kc={"zero_set_excluded":k["verified"],"restricted_action":R}
 ratio=combine(witness.dependency_token,interval_joint,kc,supply_certificate)
 return {"verified":ratio.verified,"ratio":ratio,"kernel":k,"nonlinear_premises":True}
