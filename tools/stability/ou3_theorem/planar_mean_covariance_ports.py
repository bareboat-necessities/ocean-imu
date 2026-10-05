"""Analytical Lipschitz bounds for the remaining MEKF mean<->P loop.

These are local same-history inequalities used by the planar joint-cell proof.
They deliberately do not declare a radius or forward-invariant cell.
"""
from __future__ import annotations
import math

def certificate(g=9.80665,aw_max=.02*(math.pi/10)**2,mag_norm=75.0,innovation_floor=.040976330637931824):
    # ||[x]x-[y]x||_2 <= ||x-y||.  For f=R(q)(aw-g), first-order row
    # variation is bounded by |aw-g|*attitude + aw-state variation.
    force=g+aw_max
    # Magnetic predicted vector has fixed norm 75; row variation <=75*dtheta.
    return {
      "qualification":"OU3_PLANAR_MEAN_COVARIANCE_PORT_BOUNDS_V1",
      "result_type":"PROVED analytical local inequalities; radius application OPEN",
      "acc_H_lipschitz_attitude":force,
      "acc_H_lipschitz_aw":1.0,
      "mag_H_lipschitz_attitude":mag_norm,
      "innovation_inverse_norm_upper_from_carried_floor":1.0/innovation_floor,
      "gain_identity":"K=P H^T (H P H^T+R)^-1",
      "gain_difference_bound":
        "dK=dP H^T S^-1 + P dH^T S^-1 - P H^T S^-1(dH P H^T + H dP H^T + H P dH^T)S^-1",
      "same_history_required":True,
      "joint_cell_forward_invariant":False,
      "all_time_magnetic_service_verified":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
