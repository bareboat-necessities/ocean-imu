#!/usr/bin/env python3
"""Scoped P3 execution admission consumed by regional P4.

This does not turn conditional P3 into deployment P3. It verifies that the
canonical COMPLETE_BRMM_NORMAL_LIVE_WORD itself carries every execution
premise required by the already-closed conditional P3 matrix implication.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_riccati_metric_p3 as P3


def build():
    p=P3.build()
    vf=P3.validate(p)
    if vf: raise RuntimeError('P3 validation failed: '+repr(vf))
    e=p['explicit_execution_premises']
    mandatory=e.get('mandatory_preconditions',{}) if isinstance(e,dict) else {}
    mandatory_closed=all(v is True for v in mandatory.values()) if mandatory else True
    scoped=bool(
      p['P3_CONDITIONAL_BRMM_PASS'] and p['P3_FULL_MATRIX_COMPARISON_CLOSED']
      and p['P3_FULL_WORD_ENCLOSED'] and p['actual_applied_per_axis_RS_consumed']
      and p['all_due_S_updates_required'] and p['all_valid_accelerometer_updates_required']
      and p['closed_projection_covariance_comparison_covered']
      and p['same_primitive_root_drives_entire_execution'] and mandatory_closed
      and float(p['useful_gate'])==1e-18)
    return {
      'qualification':'OU3_P4_SCOPED_CONDITIONAL_P3_EXECUTION_ADMISSION_V1',
      'canonical_source':p['canonical_source'],'P3_delta':p['useful_gate'],
      'conditional_P3_pass':p['P3_CONDITIONAL_BRMM_PASS'],
      'full_word_enclosed':p['P3_FULL_WORD_ENCLOSED'],
      'full_matrix_comparison_closed':p['P3_FULL_MATRIX_COMPARISON_CLOSED'],
      'actual_applied_per_axis_RS_consumed':p['actual_applied_per_axis_RS_consumed'],
      'all_due_S_updates_required':p['all_due_S_updates_required'],
      'all_valid_accelerometer_updates_required':p['all_valid_accelerometer_updates_required'],
      'projection_covariance_identity_covered':p['closed_projection_covariance_comparison_covered'],
      'same_primitive_root_drives_entire_execution':p['same_primitive_root_drives_entire_execution'],
      'mandatory_preconditions_closed':mandatory_closed,
      'P3_SCOPED_EXECUTION_ADMISSION_PASS':scoped,
      'P3_DEPLOYMENT_PASS':p['P3_DEPLOYMENT_PASS'],
      'global_physical_deployment_left_inclusion_closed':False,
      'P4_may_consume_conditional_P3':scoped,
    }


def validate(d):
    f=[]
    for k in ('conditional_P3_pass','full_word_enclosed','full_matrix_comparison_closed','actual_applied_per_axis_RS_consumed','all_due_S_updates_required','all_valid_accelerometer_updates_required','projection_covariance_identity_covered','same_primitive_root_drives_entire_execution','mandatory_preconditions_closed','P3_SCOPED_EXECUTION_ADMISSION_PASS','P4_may_consume_conditional_P3'):
        if d.get(k) is not True:f.append(k+' not true')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    if d.get('P3_DEPLOYMENT_PASS') is not False or d.get('global_physical_deployment_left_inclusion_closed') is not False:f.append('deployment falsely admitted')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
