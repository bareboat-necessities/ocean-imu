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


def ag_trial_action(*,T=1.0,gyro_density=.00135,bg_rw=1e-10,
                    acc_count=251,acc_h=30.0,acc_R=.04,
                    mag_count=26,mag_h=75.0,mag_R=.64):
    """Explicit linear endpoint path theta=t/T*theta1, bg=t/T*bg1.

    Dynamics theta_dot=-bg+u_g, bg_dot=u_b. Choosing this affine path gives
    u_b=bg1/T and u_g=theta1/T+(t/T)bg1. Integrating the process action and
    adding maximal direct measurement information yields an upper endpoint
    action. It is deliberately more informative than literal rows; therefore
    it is usable only after the Schur/oracle comparison is discharged.
    """
    import numpy as np
    sg2=gyro_density**2; sb2=bg_rw**2
    A=np.zeros((2,2))
    # integral (theta/T + t bg/T)^2/sg2
    A[0,0]+=1/(T*sg2);A[0,1]+=1/(2*sg2);A[1,0]+=1/(2*sg2);A[1,1]+=T/(3*sg2)
    A[1,1]+=1/(T*sb2)
    # theta(t_i)=(t_i/T) theta1; worst sum <= count for endpoint coefficient.
    A[0,0]+=acc_count*acc_h**2/acc_R+mag_count*mag_h**2/mag_R
    return variational_lower_from_action(A)

def certificate():
    counts=event_counts(); ag=ag_trial_action(acc_count=counts["acc_max"],mag_count=counts["mag_max"])
    return {
      "qualification":"OU3_PLANAR_SCHEDULER_PHASE_LOWER_CELL_V1",
      "representation":"one-second variational endpoint action, separately on exact 12/9 parity blocks",
      "scheduler_phase_uniform_event_counts":counts,
      "scheduler_phase_handling":"use maximal S-event information in the action upper bound; this is valid for every initial elapsed in [0,T_S)",
      "correction_information_handling":"literal positive-noise measurement action; accel nuisance retained by complete parity-block path, magnetic rows retained separately when constructing the covariance floor",
      "process_handling":"literal AG, LIN-OU and BA process action over the same one-second path; no one-step scalar process floor substituted",
      "ag_trial_covariance_floor":ag,
      "ag_trial_role":"conditional upper-action feasibility pending literal Schur/oracle dominance; not promoted as lower cell",
      "even_action_matrix_verified":False,
      "odd_action_matrix_verified":False,
      "lower_cell_verified":False,
      "reason":"planar parity-block comparison-path action matrices not yet evaluated",
      "theorem_closed":False,
    }

if __name__=="__main__":print(json.dumps(certificate(),indent=2,sort_keys=True))
