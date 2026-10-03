"""Non-promoting candidate for the planar scheduler-phase service enclosure.

This is the required feasibility stage before interval promotion. It uses the
fresh literal replay only to choose a candidate cell and an intentionally more
informative direct-attitude oracle to measure available magnetic-service margin.
It returns verified=False until a same-history interval map proves self-inclusion.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
ART=ROOT/"reports/results/ou3_stability/moving-compatibility-carried-1200s.json"

def oracle_service_floor(*,h_acc=30.0,p_theta=5e-4,p_bg=1.1e-6,
                         dt=.005,B=75.0,Rmag=.8**2,Racc=.2**2,
                         gyro_density=.00135,bg_rw=1e-10):
    """2x2 heading/BG conditional-information feasibility calculation.

    Direct isotropic attitude observations at h_acc/sqrt(Racc) are deliberately
    MORE informative than the actual accelerometer once nuisance is known.  A
    rigorous proof still has to establish the oracle comparison and h_acc bound.
    """
    P=np.diag([p_theta,p_bg]); X=np.array([[1.,0.],[0.,.02]]); I=np.zeros((2,2))
    qg=gyro_density**2
    for k in range(200):
        F=np.array([[1.,dt],[0.,1.]])
        Q=np.array([[qg*dt+bg_rw*dt**3/3,bg_rw*dt**2/2],
                    [bg_rw*dt**2/2,bg_rw*dt]])
        P=F@P@F.T+Q;X=F@X
        H=np.array([[h_acc,0.]])
        S=float((H@P@H.T)[0,0]+Racc);K=(P@H.T/S).reshape(2)
        A=np.eye(2)-K[:,None]@H;X=A@X;P=A@P@A.T+np.outer(K,K)*Racc
        if (k+1)%8==0:
            Hm=np.array([[B,0.]])
            S=float((Hm@P@Hm.T)[0,0]+Rmag); y=(Hm@X).reshape(2)
            I+=np.outer(y,y)/S
            K=(P@Hm.T/S).reshape(2);A=np.eye(2)-K[:,None]@Hm
            X=A@X;P=A@P@A.T+np.outer(K,K)*Rmag
    return float(np.linalg.eigvalsh(I)[0])

def scheduler_phase_invariant(elapsed,period):\n    """Literal structural invariant: valid scheduler state remains 0<=elapsed<period."""\n    return period>0 and elapsed>=0 and elapsed<period\n\ndef candidate():
    a=json.loads(ART.read_text());n=a["native"]
    floor=oracle_service_floor()
    return {
      "qualification":"OU3_PLANAR_SCHEDULER_PHASE_CELL_CANDIDATE_V1",
      "result_type":"finite feasibility candidate; NOT interval theorem",
      "candidate_cell":{
        "parity_even_spectral_upper":.04,
        "parity_odd_spectral_upper":.035,
        "attitude_covariance_upper":1.5e-5,
        "gyro_bias_block_norm_upper":5e-8,
        "attitude_bg_cross_norm_upper":1.5e-7,
        "predicted_specific_force_norm_upper":12.0,
        "scheduler_phase_interval":[0.0,float(n["period"])],
      },
      "fresh_tail_observations":{
        "parity_even_eig_max":n["parity_even_eig_max"],
        "parity_odd_eig_max":n["parity_odd_eig_max"],
        "parity_off_fro_max":n["parity_off_fro_max"],
        "attitude_covariance_eig_max":n["tail_pth_max"],
        "gyro_bias_block_norm_max":n["tail_pbg_norm_max"],
        "attitude_bg_cross_norm_max":n["tail_pth_bg_norm_max"],
        "predicted_specific_force_norm_max":n["tail_fhat_max"],
        "all_sample_root_tail_service_min":n["sliding_service_min"],
        "sample_root_windows":n["sliding_service_windows"],
        "twenty_second_full_P_max_abs_drift":n["cycle_P_max_abs_diff"],
        "twenty_second_scheduler_elapsed_abs_drift":n["cycle_scheduler_elapsed_abs_diff"],
      },
      "oracle_feasibility":{
        "h_acc":30.0,"P_theta":5e-4,"P_bg":1.1e-6,
        "magnetic_information_floor":floor,"required_floor":1.0,
        "margin":floor-1.0,
      },
      "scheduler_phase_coordinate_invariant":True,
      "verified":False,
      "open_dependencies":[
        "outward same-history propagation of the 12-state and 9-state parity covariance cells",
        "causal scheduler-phase branch enclosure over [0,T_S]",
        "all-time mean/predicted-force bound or equivalent direct bound on accelerometer attitude information",
        "formal conditional-information dominance of the nuisance-known accelerometer oracle",
        "interval propagation of each literal +/- 2x2 service Gramian at every placed root",
      ],
      "failure_class_if_cell_does_not_close":"D",
      "theorem_closed":False,
    }

if __name__=="__main__": print(json.dumps(candidate(),indent=2,sort_keys=True))
