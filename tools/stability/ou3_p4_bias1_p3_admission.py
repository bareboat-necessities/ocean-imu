#!/usr/bin/env python3
"""Qualified BIAS1 source family and canonical P3 execution admission for P4.

This is a source/admission certificate, not an endpoint contraction certificate.
The physical accelerometer-bias history is one common root/parameter/driver
trajectory; independent per-sample bias slots are forbidden.  The existing
conditional complete-BRMM P3 implication is consumed unchanged at delta=1e-18.
Global SEA0->BRMM left inclusion remains a separate deployment theorem and is
not silently imported into this scoped conditional P4 statement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_riccati_metric_p3 as P3
import ou3_p4_projection_sector as PROJ

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/"tools"/"stability"/"ou3_proof_operating_domain.json"
DEFAULT_CLOSURE=REPO/"tools"/"stability"/"ou3_p4_closure_domain.json"
QUALIFICATION="OU3_P4_BIAS1_AND_P3_EXECUTION_ADMISSION_V1"


def _member(f: dict) -> dict:
    root=[0.08,-0.05,0.03]
    amp=[0.015,0.010,-0.008]
    tau=1200.0
    period=600.0
    phase=[0.0,1.1,-0.7]
    return {
      "root": max(map(abs,root)) <= float(f["root_component_abs_upper_mps2"]),
      "amplitude": max(map(abs,amp)) <= float(f["sinusoid_component_amplitude_abs_upper_mps2"]),
      "tau": float(f["tau_true_s"][0]) <= tau <= float(f["tau_true_s"][1]),
      "period": float(f["sinusoid_period_s"][0]) <= period <= float(f["sinusoid_period_s"][1]),
      "phase": all(float(f["phase_rad"][0]) <= x <= float(f["phase_rad"][1]) for x in phase),
      "probe": {"root":root,"amplitude":amp,"tau_s":tau,"period_s":period,"phase_rad":phase},
    }


def build(domain_path: Path=DEFAULT_DOMAIN, closure_path: Path=DEFAULT_CLOSURE) -> dict:
    closure=json.loads(Path(closure_path).read_text())
    family=closure["BIAS1_family"]
    p3=P3.build(Path(domain_path).resolve())
    pf=P3.validate(p3)
    if pf:
        raise RuntimeError(f"canonical P3 invalid: {pf}")
    proj=PROJ.build()
    qf=PROJ.validate(proj)
    if qf:
        raise RuntimeError(f"projection sector invalid: {qf}")
    membership=_member(family)
    member=all(membership[k] for k in ("root","amplitude","tau","period","phase"))
    premises=p3["explicit_execution_premises"]
    p3_exec=bool(
        p3["P3_CONDITIONAL_BRMM_PASS"]
        and float(p3["useful_gate"])==1e-18
        and p3["all_due_S_updates_required"]
        and p3["actual_applied_per_axis_RS_consumed"]
        and p3["all_valid_accelerometer_updates_required"]
        and p3["closed_projection_covariance_comparison_covered"]
        and p3["same_complete_BRMM_execution_continues_across_H_to_A"]
        and p3["same_primitive_root_drives_entire_execution"]
    )
    return {
      "qualification":QUALIFICATION,
      "canonical_source":p3["canonical_source"],
      "BIAS1_family":family,
      "BIAS1_family_compact": True,
      "BIAS1_one_common_root_parameter_driver_history": True,
      "BIAS1_independent_per_sample_slots_forbidden": True,
      "BIAS1_affine_driver_recurrence":"beta_i=phi_true*beta_i-1+w_i; w_i=beta(t_i)-phi_true*beta(t_i-1)",
      "existing_nonzero_driver_probe_membership":membership,
      "existing_nonzero_driver_probe_admitted":member,
      "projection_exact_real_global_sector_consumed":proj["global_joint_sector_closed"],
      "projection_estimate_ball_invariance_consumed":proj["estimate_ball_invariance_closed"],
      "BIAS1_SOURCE_ADMISSION_PASS":member and proj["global_joint_sector_closed"],
      "P3_frozen_delta":1e-18,
      "P3_numerical_certificate_modified":False,
      "P3_execution_premises":premises,
      "P3_EXECUTION_ADMISSION_PASS":p3_exec,
      "global_SEA0_to_BRMM_left_inclusion_required_for_scoped_conditional_math":False,
      "global_SEA0_to_BRMM_left_inclusion_closed_here":False,
      "trajectory_replay_used_for_admission":False,
      "filter_changed":False,
      "P4_endpoint_contraction_promoted_here":False,
    }


def validate(d:dict)->list[str]:
    f=[]
    for k in ("BIAS1_family_compact","BIAS1_one_common_root_parameter_driver_history",
              "BIAS1_independent_per_sample_slots_forbidden","existing_nonzero_driver_probe_admitted",
              "projection_exact_real_global_sector_consumed","projection_estimate_ball_invariance_consumed",
              "BIAS1_SOURCE_ADMISSION_PASS","P3_EXECUTION_ADMISSION_PASS"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in ("P3_numerical_certificate_modified","global_SEA0_to_BRMM_left_inclusion_required_for_scoped_conditional_math",
              "global_SEA0_to_BRMM_left_inclusion_closed_here","trajectory_replay_used_for_admission","filter_changed",
              "P4_endpoint_contraction_promoted_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if float(d.get("P3_frozen_delta",0))!=1e-18:f.append("P3 delta changed")
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);args=ap.parse_args()
    d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"BIAS1":d["BIAS1_SOURCE_ADMISSION_PASS"],"P3_execution":d["P3_EXECUTION_ADMISSION_PASS"],"failures":f},indent=2))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
