#!/usr/bin/env python3
"""Single fail-closed conditional mathematical P4 promotion gate.

The mathematical theorem may consume explicit execution premises, just as the
canonical P3 theorem does. Target-toolchain/device qualification is reported
separately and is never inferred from mathematical P4.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1
import ou3_p4_p3_execution_admission as P3A
import ou3_p4_exact_chord_signed_master_bridge as BRIDGE
import ou3_p4_kalman_reset_binary32_iss as FP
QUALIFICATION='OU3_P4_FINAL_FAIL_CLOSED_GATE_V3'

def build():
    e=ENTRY.build();b=BIAS1.build();p=P3A.build();g=BRIDGE.build();fp=FP.build()
    bad={'entry':ENTRY.validate(e),'bias1':BIAS1.validate(b),'p3':P3A.validate(p),'bridge':BRIDGE.validate(g),'fp':FP.validate(fp)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('final P4 prerequisite validation failed: '+repr(bad))
    admissions=bool(e['P4_may_use_as_regional_entry_hypothesis'] and b['BIAS1_SOURCE_ADMISSION_PASS'] and p['P3_SCOPED_EXECUTION_ADMISSION_PASS'])
    backbone=bool(g['universal_full_entry_finite_angle_differential_backbone_closed'])
    nonlinear_graph=bool(g['exact_chord_graph_ready_for_augmented_master'] and g['joint_ISS_master_ready_for_BIAS1_and_roundoff'] and g['physical_BIAS1_projection_and_coefficient_prerequisites_ready'])
    endpoint=bool(g['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'])
    prefixes=bool(g['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'])
    hard_prefix=bool(g['same_graph_every_prefix_hard_domain_retention_closed'])
    arithmetic=bool(fp['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally'] and fp['additive_ISS_channel_complete_for_conditional_P4'])
    platform_qualified=bool(fp['deployment_finite_precision_qualification_closed'])
    motion=bool(admissions and backbone and nonlinear_graph and endpoint and prefixes and hard_prefix and arithmetic)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','declared_domain_shrunk':False,
      'filter_changed':False,'quality_gates_changed':False,'P3_delta':1e-18,
      'hard_entry_set_admitted_without_covariance_membership':admissions,
      'universal_full_entry_finite_angle_differential_backbone_closed':backbone,
      'H18_finite_angle_worst_LDLT_pivot_lower':g['universal_H18_worst_LDLT_pivot_lower'],
      'A21_first_active_ba_margin_lower':g['universal_A21_first_active_ba_margin_lower'],
      'exact_chord_projection_BIAS1_same_history_graph_ready':nonlinear_graph,
      'source_uniform_endpoint_augmented_LDLT_closed':endpoint,
      'source_uniform_every_prefix_augmented_LDLT_closed':prefixes,
      'source_uniform_every_prefix_hard_domain_retention_closed':hard_prefix,
      'conditional_full_shipping_finite_precision_additive_ISS_closed':arithmetic,
      'target_toolchain_finite_precision_qualified':platform_qualified,
      'P4_DEPLOYMENT_PASS':bool(motion and platform_qualified),
      'point_capture_can_promote':False,'rowwise_coefficient_boxes_can_promote':False,'differential_backbone_alone_can_promote':False,
      'P4_MOTION_PASS':motion,'P4_PASS':motion,'P5_MAY_START':motion,
      'remaining_mathematical_P4_blockers':[x for x,ok in (
        ('source-uniform exact-graph endpoint augmented LDLT',endpoint),
        ('source-uniform exact-graph every-prefix augmented LDLT',prefixes),
        ('same exact graph every-prefix hard-domain retention',hard_prefix)) if not ok],
      'remaining_deployment_blockers':[x for x,ok in (
        ('target-toolchain qualification of the declared binary32 LDLT/libm forward postconditions',platform_qualified),) if not ok]}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    for k in ('hard_entry_set_admitted_without_covariance_membership','universal_full_entry_finite_angle_differential_backbone_closed','exact_chord_projection_BIAS1_same_history_graph_ready','conditional_full_shipping_finite_precision_additive_ISS_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('declared_domain_shrunk','filter_changed','quality_gates_changed','source_uniform_endpoint_augmented_LDLT_closed','source_uniform_every_prefix_augmented_LDLT_closed','source_uniform_every_prefix_hard_domain_retention_closed','target_toolchain_finite_precision_qualified','P4_DEPLOYMENT_PASS','point_capture_can_promote','rowwise_coefficient_boxes_can_promote','differential_backbone_alone_can_promote','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for k in ('H18_finite_angle_worst_LDLT_pivot_lower','A21_first_active_ba_margin_lower'):
        if float(d.get(k,0))<=0:f.append(k+' not positive')
    if len(d.get('remaining_mathematical_P4_blockers',[]))!=3:f.append('expected three mathematical blockers')
    if len(d.get('remaining_deployment_blockers',[]))!=1:f.append('expected one deployment blocker')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
