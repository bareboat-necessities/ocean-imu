#!/usr/bin/env python3
"""Scoped conditional P3 execution admission consumed by regional P4.

This does not turn conditional P3 into deployment P3. It checks the actual
premise manifest bound to the declared operating-domain hash and verifies that
the canonical COMPLETE_BRMM_NORMAL_LIVE_WORD P3 certificate consumes those
premises unchanged at delta=1e-18. Runtime ``Live`` alone is explicitly not an
admission test.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_riccati_metric_p3 as P3
import ou3_brmm_p3_premises as PREMISES


def _execution_manifest_closed(e: dict) -> bool:
    x=e["execution_premises"]
    bool_keys=(
        "same_complete_physical_frontend_tuner_geometry_history",
        "declared_Normal_Live_branch_not_just_runtime_Live_flag",
        "lever_arm_disabled_and_vibration_guard_dormant_transparent",
        "accepted_accelerometer_at_every_valid_IMU_sample",
        "no_accelerometer_rejection_in_word",
        "PE_geometry_attached_to_actual_measurement_jacobians",
        "committed_OU_parameters_obey_shipping_clamp_and_commit_invariants",
        "every_due_S_with_actual_anisotropic_RS_and_full_cross_covariance",
        "shipping_configured_R_full_Q_and_PSD_aw_floors",
        "successful_Joseph_updates_and_immediate_covariance_resets",
        "no_hard_attitude_rewrite_inside_same_mode_word",
        "source_generated_Live_seed_and_shipping_H_to_A_release",
        "premises_continue_after_every_projection_and_hybrid_entry",
    )
    if not all(x.get(k) is True for k in bool_keys):
        return False
    pe=x["vector_PE"]
    positive=(
        "vector_pe_recurrence_window_s","specific_force_norm_lower_mps2",
        "specific_force_norm_upper_mps2","magnetic_vector_norm_lower_uT",
        "magnetic_vector_norm_upper_uT","vector_sine_separation_lower",
        "body_rate_norm_upper_deg_s",
    )
    if not all(math.isfinite(float(pe[k])) and float(pe[k])>0.0 for k in positive):
        return False
    if float(pe["specific_force_norm_lower_mps2"]) > float(pe["specific_force_norm_upper_mps2"]):
        return False
    if float(pe["magnetic_vector_norm_lower_uT"]) > float(pe["magnetic_vector_norm_upper_uT"]):
        return False
    return True


def build():
    p=P3.build()
    vf=P3.validate(p)
    if vf: raise RuntimeError('P3 validation failed: '+repr(vf))
    e=PREMISES.build()
    ef=PREMISES.validate(e)
    if ef: raise RuntimeError('P3 premise validation failed: '+repr(ef))
    manifest_closed=_execution_manifest_closed(e)
    scoped=bool(
      manifest_closed
      and p['P3_CONDITIONAL_BRMM_PASS'] and p['P3_FULL_MATRIX_COMPARISON_CLOSED']
      and p['P3_FULL_WORD_ENCLOSED'] and p['actual_applied_per_axis_RS_consumed']
      and p['all_due_S_updates_required'] and p['all_valid_accelerometer_updates_required']
      and p['closed_projection_covariance_comparison_covered']
      and p['same_primitive_root_drives_entire_execution']
      and p['explicit_execution_premises']['domain_sha256']==e['domain_sha256']
      and float(p['useful_gate'])==1e-18 and float(e['delta'])==1e-18)
    return {
      'qualification':'OU3_P4_SCOPED_CONDITIONAL_P3_EXECUTION_ADMISSION_V2',
      'canonical_source':p['canonical_source'],'P3_delta':p['useful_gate'],
      'premise_domain_sha256':e['domain_sha256'],
      'conditional_P3_pass':p['P3_CONDITIONAL_BRMM_PASS'],
      'full_word_enclosed':p['P3_FULL_WORD_ENCLOSED'],
      'full_matrix_comparison_closed':p['P3_FULL_MATRIX_COMPARISON_CLOSED'],
      'actual_applied_per_axis_RS_consumed':p['actual_applied_per_axis_RS_consumed'],
      'all_due_S_updates_required':p['all_due_S_updates_required'],
      'all_valid_accelerometer_updates_required':p['all_valid_accelerometer_updates_required'],
      'projection_covariance_identity_covered':p['closed_projection_covariance_comparison_covered'],
      'same_primitive_root_drives_entire_execution':p['same_primitive_root_drives_entire_execution'],
      'execution_manifest_closed':manifest_closed,
      'runtime_Live_flag_implies_all_premises':e['runtime_Live_flag_implies_all_premises'],
      'BRMM_alone_implies_vector_PE':e['BRMM_alone_implies_vector_PE'],
      'physical_execution_admission_proved_here':e['physical_execution_admission_proved_here'],
      'P3_SCOPED_EXECUTION_ADMISSION_PASS':scoped,
      'P3_DEPLOYMENT_PASS':p['P3_DEPLOYMENT_PASS'],
      'global_physical_deployment_left_inclusion_closed':False,
      'P4_may_consume_conditional_P3':scoped,
    }


def validate(d):
    f=[]
    for k in ('conditional_P3_pass','full_word_enclosed','full_matrix_comparison_closed',
              'actual_applied_per_axis_RS_consumed','all_due_S_updates_required',
              'all_valid_accelerometer_updates_required','projection_covariance_identity_covered',
              'same_primitive_root_drives_entire_execution','execution_manifest_closed',
              'P3_SCOPED_EXECUTION_ADMISSION_PASS','P4_may_consume_conditional_P3'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('runtime_Live_flag_implies_all_premises','BRMM_alone_implies_vector_PE',
              'physical_execution_admission_proved_here','P3_DEPLOYMENT_PASS',
              'global_physical_deployment_left_inclusion_closed'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
