"""Scheduler-phase-uniform lower covariance cell: structural certificate.

The exact lower cell is represented variationally rather than as entrywise P
intervals. For a one-second horizon, posterior precision at the endpoint is the
minimum action of process controls plus all literal measurements. Any explicit
comparison path gives an UPPER action matrix A, hence P_end >= A^-1.

This module assembles the proof obligations for the planar 12+9 parity blocks
and proves the scheduler-phase event-count part. Numeric action matrices remain
fail-closed until the planar Hermite/control action is evaluated.
"""
from __future__ import annotations
import math,json

EVEN_DIM,ODD_DIM=12,9

def event_counts(window_s=1.0,dt_min=.004,mag_period=.04,Ts_min=.12,Ts_max=.18):
    # Acc every IMU sample. For a placed one-second window with arbitrary
    # scheduler phase, the largest number of due S events is ceil(T/Ts_min)+1.
    return {"acc_max":math.ceil(window_s/dt_min)+1,
            "mag_min":math.floor(window_s/mag_period),
            "mag_max":math.ceil(window_s/mag_period)+1,
            "S_max":math.ceil(window_s/Ts_min)+1,
            "S_min":max(0,math.floor(window_s/Ts_max)-1)}

def variational_lower_from_action(action):
    """Return smallest covariance floor from SPD endpoint action."""
    import numpy as np
    A=np.asarray(action,float); ev=np.linalg.eigvalsh((A+A.T)/2)
    if ev[0]<=0: raise ValueError("action must be SPD")
    return {"action_lambda_max":float(ev[-1]),"covariance_scalar_floor":float(1/ev[-1]),
            "covariance_matrix":np.linalg.inv(A).tolist()}

def certificate():
    counts=event_counts()
    return {
      "qualification":"OU3_PLANAR_SCHEDULER_PHASE_LOWER_CELL_V1",
      "representation":"one-second variational endpoint action, separately on exact 12/9 parity blocks",
      "scheduler_phase_uniform_event_counts":counts,
      "scheduler_phase_handling":"use maximal S-event information in the action upper bound; this is valid for every initial elapsed in [0,T_S)",
      "correction_information_handling":"literal positive-noise measurement action; accel nuisance retained by complete parity-block path, magnetic rows retained separately when constructing the covariance floor",
      "process_handling":"literal AG, LIN-OU and BA process action over the same one-second path; no one-step scalar process floor substituted",
      "even_action_matrix_verified":False,
      "odd_action_matrix_verified":False,
      "lower_cell_verified":False,
      "reason":"planar parity-block comparison-path action matrices not yet evaluated",
      "theorem_closed":False,
    }

if __name__=="__main__":print(json.dumps(certificate(),indent=2,sort_keys=True))
