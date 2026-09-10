#!/usr/bin/env python3
"""Bind non-S same-cell reset domains and defer S=0 to the centered prefix graph.

Accelerometer/magnetometer correction ceilings come from the exact same-cell
Joseph covariance identity, never a rowwise K box.  The S=0 residual contains a
shared integration-origin cancellation and therefore cannot be bounded from the
legacy e_S hard radius.  Its reset-chart sector is a production augmented-prefix
obligation using d=K(e_S-S_true).  This module remains fail-closed until that
integrated target closes.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_p4_same_cell_correction_domain as CORR
import ou3_p4_reset_signed_absorption_diagnostic as SECTOR
import ou3_p4_reset_graph_iqc as RIQC
import ou3_p4_exact_reset_joint_coordinate as RJOIN
import ou3_p4_exact_reset_lifted_qc as RLIFT

QUALIFICATION='OU3_P4_SAME_CELL_RESET_DOMAIN_BINDING_V3'

def build():
    corr=CORR.build();sector0=SECTOR.build();riqc=RIQC.build();rjoin=RJOIN.build()
    bad={'correction':CORR.validate(corr),'sector':SECTOR.validate(sector0),'reset_iqc':RIQC.validate(riqc),'reset_joint':RJOIN.validate(rjoin)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset-domain binding prerequisites failed: '+repr(bad))
    modes={};vector_closed=True
    for mode,m in corr['modes'].items():
        delta=float(m['bounded_vector_event_correction_norm_upper']);inside=bool(m['vector_events_reset_utility_domain_closed'])
        lift=None;lf=[];sec=None;mu=None;valid=False
        if inside:
            sec=SECTOR.homogeneous_sector(float(corr['residual_bounds']['attitude_cayley_norm_upper']),delta);mu=float(sec['reset_defect_over_correction_norm_upper'])
            lift=RLIFT.build_for_same_cell(delta);lf=RLIFT.validate(lift)
            valid=bool(sec['chart_safe'] and sec['homogeneous_sector_dominates_endpoint_absolute_bound'] and not lf and lift['same_cell_correction_domain_supplied'] and not lift['rowwise_K_correction_domain_used'])
        vector_closed=vector_closed and inside and valid
        modes[mode]={'bounded_vector_event_correction_norm_upper':delta,'limiting_bounded_vector_event':m['limiting_bounded_vector_event'],
          'all_event_bounds':m['events'],'bounded_vector_events_inside_exact_reset_utility_domain':inside,
          'bounded_vector_event_homogeneous_reset_sector':sec,'bounded_vector_event_reset_defect_over_correction_norm_upper':mu,
          'bounded_vector_event_same_cell_exact_reset_lift':lift,'bounded_vector_event_same_cell_exact_reset_lift_failures':lf,
          'bounded_vector_event_uniform_reset_gain':valid,'S_zero_centered_prefix_graph_required':True,
          'S_zero_uniform_reset_gain_closed_here':False,'uniform_reset_gain_over_all_events':False,
          # Compatibility fields consumed by the final bridge/reporting.  They
          # intentionally remain unresolved because S=0 is now integrated.
          'correction_norm_upper':None,'limiting_event':'S_zero:centered-prefix-required',
          'inside_exact_reset_utility_domain':False,'homogeneous_reset_sector':None,
          'reset_defect_over_correction_norm_upper':None,'same_cell_exact_reset_lift':None,
          'same_cell_exact_reset_lift_failures':[],'uniform_reset_gain_over_zero_to_ceiling':False}
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_Joseph_correction_domain_consumed':True,'rowwise_K_bound_used':False,'independent_K_box_used':False,
      'exact_reset_joint_coordinate_consumed':True,'same_cell_parameterized_exact_reset_lift_builder_used':True,'same_cell_d_equals_Etheta_K_y_retained':True,
      'parameterized_dense_reset_IQC_consumed':True,'parameterized_dense_reset_IQC_consumed_for_bounded_vector_events':True,
      'scalar_delta_used_only_for_reset_chart_sector_validity_and_same_graph_QC_compactness':True,
      'scalar_delta_used_as_independent_storage_port':False,'historical_rowwise_reset_lift_path_consumed':False,
      'legacy_eS_radius_used_for_S_zero':False,'S_zero_reset_chart_must_close_inside_centered_augmented_prefix_graph':True,
      'bounded_vector_event_reset_domain_closed':vector_closed,'modes':modes,
      'source_uniform_same_graph_correction_domain_closed':False,'exact_reset_lifted_QCs_consumed_from_same_cell_domain':False,
      'homogeneous_same_cell_reset_IQC_consumed':False,'reset_gain_uniform_over_zero_to_correction_ceiling':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'every_prefix_hard_domain_retention_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_Joseph_correction_domain_consumed','exact_reset_joint_coordinate_consumed','same_cell_parameterized_exact_reset_lift_builder_used','same_cell_d_equals_Etheta_K_y_retained','parameterized_dense_reset_IQC_consumed','parameterized_dense_reset_IQC_consumed_for_bounded_vector_events','scalar_delta_used_only_for_reset_chart_sector_validity_and_same_graph_QC_compactness','S_zero_reset_chart_must_close_inside_centered_augmented_prefix_graph'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_bound_used','independent_K_box_used','scalar_delta_used_as_independent_storage_port','historical_rowwise_reset_lift_path_consumed','legacy_eS_radius_used_for_S_zero','source_uniform_same_graph_correction_domain_closed','exact_reset_lifted_QCs_consumed_from_same_cell_domain','homogeneous_same_cell_reset_IQC_consumed','reset_gain_uniform_over_zero_to_correction_ceiling','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','every_prefix_hard_domain_retention_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,m in d.get('modes',{}).items():
        x=float(m.get('bounded_vector_event_correction_norm_upper',math.nan))
        if not(math.isfinite(x) and x>=0):f.append(mode+' vector correction ceiling invalid')
        if m.get('S_zero_centered_prefix_graph_required') is not True or m.get('S_zero_uniform_reset_gain_closed_here') is not False:f.append(mode+' S-zero scope invalid')
        if m.get('correction_norm_upper') is not None:f.append(mode+' all-event scalar ceiling reintroduced')
    return f

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'vector_reset_domain':d['bounded_vector_event_reset_domain_closed'],'S_zero':'integrated-prefix-required','all_events_closed':d['source_uniform_same_graph_correction_domain_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
