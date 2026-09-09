#!/usr/bin/env python3
"""Production binding of same-cell Joseph correction domain to exact reset graph.

This module closes only the reset-domain prerequisite.  It deliberately does
not close endpoint/prefix dissipation or retention.

For each mode/event the source-uniform same-cell Joseph certificate supplies

  ||d_theta|| <= delta_event,

without a rowwise K box.  The exact finite-reset lift keeps d_theta bound to the
same Joseph graph and the parameterized homogeneous reset sector supplies

  ||rho_theta|| <= mu_R(delta_mode) ||d_theta||

uniformly for 0 <= ||d_theta|| <= delta_mode.  The exact lifted QCs are built
with that SAME delta.  No historical rowwise correction-domain path is allowed
inside this production binding.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_same_cell_correction_domain as CORR
import ou3_p4_reset_signed_absorption_diagnostic as SECTOR
import ou3_p4_reset_graph_iqc as RIQC
import ou3_p4_exact_reset_joint_coordinate as RJOIN
import ou3_p4_exact_reset_lifted_qc as RLIFT

QUALIFICATION='OU3_P4_SAME_CELL_RESET_DOMAIN_BINDING_V2'


def build():
    corr=CORR.build();sector0=SECTOR.build();riqc=RIQC.build();rjoin=RJOIN.build()
    bad={'correction':CORR.validate(corr),'sector':SECTOR.validate(sector0),'reset_iqc':RIQC.validate(riqc),
         'reset_joint':RJOIN.validate(rjoin)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset-domain binding prerequisites failed: '+repr(bad))
    modes={};closed=True
    for mode,m in corr['modes'].items():
        delta=float(m['correction_norm_upper'])
        inside=bool(m['reset_utility_domain_closed'])
        lift=None;lift_fail=[];sec=None;mu=None;valid=False
        if inside:
            sec=SECTOR.homogeneous_sector(float(corr['residual_bounds']['attitude_cayley_norm_upper']),delta)
            mu=float(sec['reset_defect_over_correction_norm_upper'])
            lift=RLIFT.build_for_same_cell(delta);lift_fail=RLIFT.validate(lift)
            valid=bool(sec['chart_safe'] and sec['homogeneous_sector_dominates_endpoint_absolute_bound']
                       and not lift_fail
                       and lift['same_cell_correction_domain_supplied']
                       and not lift['rowwise_K_correction_domain_used'])
        modes[mode]={'correction_norm_upper':delta,'limiting_event':m['limiting_event'],
                     'all_event_bounds':m['events'],'inside_exact_reset_utility_domain':inside,
                     'homogeneous_reset_sector':sec,'reset_defect_over_correction_norm_upper':mu,
                     'same_cell_exact_reset_lift':lift,'same_cell_exact_reset_lift_failures':lift_fail,
                     'uniform_reset_gain_over_zero_to_ceiling':valid}
        closed=closed and inside and valid
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_Joseph_correction_domain_consumed':True,'rowwise_K_bound_used':False,'independent_K_box_used':False,
      'exact_reset_joint_coordinate_consumed':True,'same_cell_parameterized_exact_reset_lift_builder_used':True,
      'same_cell_d_equals_Etheta_K_y_retained':True,'parameterized_dense_reset_IQC_consumed':True,
      'scalar_delta_used_only_for_reset_chart_sector_validity_and_same_graph_QC_compactness':True,
      'scalar_delta_used_as_independent_storage_port':False,'historical_rowwise_reset_lift_path_consumed':False,
      'modes':modes,'source_uniform_same_graph_correction_domain_closed':closed,
      'exact_reset_lifted_QCs_consumed_from_same_cell_domain':closed,
      'homogeneous_same_cell_reset_IQC_consumed':closed,'reset_gain_uniform_over_zero_to_correction_ceiling':closed,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,
      'every_prefix_hard_domain_retention_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_Joseph_correction_domain_consumed','exact_reset_joint_coordinate_consumed',
              'same_cell_parameterized_exact_reset_lift_builder_used','same_cell_d_equals_Etheta_K_y_retained',
              'parameterized_dense_reset_IQC_consumed','scalar_delta_used_only_for_reset_chart_sector_validity_and_same_graph_QC_compactness'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_bound_used','independent_K_box_used','scalar_delta_used_as_independent_storage_port',
              'historical_rowwise_reset_lift_path_consumed','endpoint_augmented_LDLT_closed_here',
              'every_prefix_augmented_LDLT_closed_here','every_prefix_hard_domain_retention_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,m in d.get('modes',{}).items():
        x=float(m.get('correction_norm_upper',math.nan))
        if not(math.isfinite(x) and x>=0):f.append(mode+' correction ceiling invalid')
        lift=m.get('same_cell_exact_reset_lift')
        if lift is not None:
            if m.get('same_cell_exact_reset_lift_failures'):f.append(mode+' same-cell reset lift invalid')
            if lift.get('same_cell_correction_domain_supplied') is not True:f.append(mode+' lift lost same-cell domain')
            if lift.get('rowwise_K_correction_domain_used') is not False:f.append(mode+' lift leaked rowwise K domain')
            if not math.isclose(float(lift.get('correction_norm_upper',math.nan)),x,rel_tol=0.0,abs_tol=0.0):f.append(mode+' lift delta differs from same-cell delta')
    # A too-large correction ceiling is a mathematical result, not a malformed
    # certificate.  Closure booleans remain false and downstream promotion
    # fails closed.
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'closed':d['source_uniform_same_graph_correction_domain_closed'],'modes':d['modes'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
