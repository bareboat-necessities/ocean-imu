#!/usr/bin/env python3
"""All-family A21 same-graph correction/reset chart cover.

The original radial event AD cover used the BIAS1 true-bias envelope when
checking A21 correction/reset charts.  The all-family bias theorem proves that
BIAS0/BIAS1/BIAS2 share the same declared component envelope (up to outward
rounding), but the three physical driver families remain distinct.

This module makes that distinction executable: for every retained estimator
image and every exact radial cell, it reruns S=0/accelerometer/magnetometer
same-graph chart certification separately for BIAS0, BIAS1 and BIAS2 using that
family's own authoritative component envelope.  P/H/R, actual R_S, nonlinear
sectors and reset binding are inherited from the exact same event machinery.
No family is inferred from BIAS1.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_brmm_kernel_event_capture_bridge as BRIDGE
import ou3_p4_brmm_kernel_radial_event_ad_cover as BASE
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_bias_family_joint_iss_supply as BIAS

SCHEMA=1
QUALIFICATION='OU3_P4_ALL_BIAS_SAME_GRAPH_CORRECTION_RESET_CHART_COVER_V1'
P3_DELTA=1e-18

def family_box(family,bias):
    if family not in BIAS.REQUIRED_BIAS_FAMILIES:raise ValueError('unknown bias family')
    r=math.nextafter(float(bias['true_bias_component_envelope_mps2'][family]),math.inf)
    return [Interval(-r,r) for _ in range(3)]

def build():
    base=BASE.build();bf=BASE.validate(base);bias=BIAS.build();bif=BIAS.validate(bias)
    if bf or bif:raise RuntimeError(f'all-bias chart prerequisites failed base={bf} bias={bif}')
    if not bias['families_share_one_true_bias_envelope']:
        raise RuntimeError('shared component envelope lemma is not closed')
    sample,_,meta,_,images=BRIDGE.capture_smoke();parts=RADIAL.partition(2);records=[]
    for family in BIAS.REQUIRED_BIAS_FAMILIES:
      tb=family_box(family,bias)
      for image in images:
        by=BASE._captured_by_kind(meta['A_event_cells'])
        for p in parts:
          state=RADIAL.radial_state_box('A',p.interval)
          rows=[]
          for ordinal,kind in enumerate(BASE.JOSEPH,2):
            rows.append({'kind':kind,**BASE.inspect_event(image,'A',kind,state,p.interval,by[kind],sample,ordinal,tb)})
          records.append({'family':family,'image_token':image.source_token,'radial':p.interval.as_list(),'events':rows})
    allrows=[e for r in records for e in r['events']]
    closed=bool(allrows) and all(e.get('closed') is True for e in allrows)
    same=all(e.get('same_P_H_R_cell') is True for e in allrows)
    charts=all(e.get('same_graph_correction_target_closed') is True and e.get('first_exit_reset_structure_closed') is True for e in allrows)
    actual_rs=all(e.get('R_provenance')==BASE.EVENTS.ACTUAL_RS_PROVENANCE for e in allrows if e['kind']=='S_zero')
    families=set(r['family'] for r in records)==set(BIAS.REQUIRED_BIAS_FAMILIES)
    deltas=[float(e['certified_correction_delta']) for e in allrows if math.isfinite(float(e.get('certified_correction_delta',math.nan)))]
    pivots=[float(e['same_graph_chart_certificate']['minimum_pivot_lower']) for e in allrows if e.get('same_graph_chart_certificate') and math.isfinite(float(e['same_graph_chart_certificate'].get('minimum_pivot_lower',math.nan)))]
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'BIAS0_BIAS1_BIAS2_certified_separately':families,'shared_component_envelope_lemma_consumed':True,'family_driver_recurrences_collapsed_into_one':False,'BIAS0_or_BIAS2_inferred_from_BIAS1':False,'same_estimator_owned_P_H_R_retained':same,'actual_applied_RS_retained':actual_rs,'same_graph_correction_and_reset_chart_closed_all_families':bool(closed and charts),'record_count':len(records),'event_count':len(allrows),'worst_certified_correction_delta':max(deltas,default=math.nan),'minimum_chart_LDLT_pivot_lower':min(pivots,default=math.nan),'records':records,
      'source_uniform_family_driver_endpoint_LDLT_closed_here':False,'source_uniform_every_prefix_LDLT_closed_here':False,'first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'combine this all-family correction/reset chart cover with each family-specific sequential source/bias/finite-precision endpoint and every-prefix LDLT over the universal selector relation'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('BIAS0_BIAS1_BIAS2_certified_separately','shared_component_envelope_lemma_consumed','same_estimator_owned_P_H_R_retained','actual_applied_RS_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    # Chart closure is a mathematical result rather than a schema axiom: retain
    # failing evidence if a family/radial cell does not certify.
    for k in ('family_driver_recurrences_collapsed_into_one','BIAS0_or_BIAS2_inferred_from_BIAS1','source_uniform_family_driver_endpoint_LDLT_closed_here','source_uniform_every_prefix_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if int(d.get('record_count',0))<=0 or int(d.get('event_count',0))!=3*int(d.get('record_count',0)):f.append('all-family chart event count mismatch')
    if d.get('same_graph_correction_and_reset_chart_closed_all_families'):
        if not(0<float(d.get('worst_certified_correction_delta',0))<3):f.append('invalid correction delta')
        if not float(d.get('minimum_chart_LDLT_pivot_lower',0))>0:f.append('invalid chart pivot')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'records':d['record_count'],'all_family_chart':d['same_graph_correction_and_reset_chart_closed_all_families'],'worst_delta':d['worst_certified_correction_delta'],'min_pivot':d['minimum_chart_LDLT_pivot_lower'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
