#!/usr/bin/env python3
"""Continuous BRMM/radial event-lineage cover for OU-III P4.

This is the theorem-facing lineage object between the same-signal estimator and
the retained signed-master/ISS backends.

A lineage retains simultaneously:
  * one BRMM/private-Mahony estimator ancestry;
  * literal Kalman/Joseph/reset event ancestry owned by each estimator image;
  * the hard-entry radial source coordinate as CLOSED INTERVALS whose union is
    exactly [0,1], never sampled radial points;
  * reachable covariance continuity between literal events;
  * the requirement that production augmented masters consume the retained
    exact-chord signed ledger, physical BIAS1 supply and binary32 ISS channels.

This module deliberately does not fabricate a master matrix.  It fails closed
until the real same-history event master maps are supplied for every retained
lineage/prefix.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_source_cover_transition as COV
import ou3_p4_exact_chord_signed_master_bridge as SIGNED_MASTER
import ou3_p4_complete_brmm_joint_sector_master as MASTER
import ou3_p4_bias1_joint_iss_supply as BIASISS
import ou3_p4_kalman_reset_binary32_iss as FPISS
import ou3_p4_joint_brmm_prefix_ldlt as PREFIX
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_frontend_state_step as FRONT

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_CONTINUOUS_RADIAL_EVENT_LINEAGE_COVER_V1'
ZERO=Interval.point(0.0);ONE=Interval.point(1.0)

def I(x):return Interval.point(float(x))

@dataclass(frozen=True)
class RadialCell:
    token:str
    interval:Interval
    parent_token:str|None

@dataclass(frozen=True)
class EventLineage:
    radial:RadialCell
    estimator_token:str
    cells:tuple[COVER.SourceCoverCell,...]


def root_radial(token='RADIAL_ROOT'):
    return RadialCell(token,Interval(0.0,1.0),None)

def split_radial(cell:RadialCell,*,split:float|None=None):
    lo,hi=cell.interval.lo,cell.interval.hi
    if not (0.0<=lo<hi<=1.0):raise ValueError('radial cell must be nondegenerate subset of [0,1]')
    m=float(split) if split is not None else lo+(hi-lo)*0.5
    if not(lo<m<hi):raise ValueError('split must be strictly interior')
    # Shared exact binary64 boundary means union(children) is exactly parent;
    # no uncovered numerical crack and no radial point sampling is introduced.
    return (RadialCell(cell.token+':L',Interval(lo,m),cell.token),RadialCell(cell.token+':R',Interval(m,hi),cell.token))

def partition_covers_unit(cells:Sequence[RadialCell]):
    if not cells:return False
    ordered=sorted(cells,key=lambda c:(c.interval.lo,c.interval.hi))
    if ordered[0].interval.lo!=0.0 or ordered[-1].interval.hi!=1.0:return False
    if any(c.interval.lo<0 or c.interval.hi>1 or c.interval.lo>c.interval.hi for c in ordered):return False
    return all(ordered[i].interval.hi==ordered[i+1].interval.lo for i in range(len(ordered)-1))

def refine_partition(cells:Sequence[RadialCell],token:str):
    out=[];found=False
    for c in cells:
        if c.token==token:
            out.extend(split_radial(c));found=True
        else:out.append(c)
    if not found:raise KeyError(token)
    if not partition_covers_unit(out):raise RuntimeError('radial refinement lost exact [0,1] cover')
    return tuple(out)

def validate_event_lineage(lineage:EventLineage):
    f=[];cells=lineage.cells
    if not cells:return ['empty event lineage']
    if lineage.radial.interval.lo<0 or lineage.radial.interval.hi>1:f.append('radial cell outside [0,1]')
    for i,c in enumerate(cells):
        vf=COVER.validate_cell(c,require_estimator_provenance=True)
        f.extend(f'event {i}: {x}' for x in vf)
        if c.estimator_source_token!=lineage.estimator_token:f.append(f'event {i}: estimator token detached')
        if c.radial_scale is None or c.radial_scale.lo!=lineage.radial.interval.lo or c.radial_scale.hi!=lineage.radial.interval.hi:f.append(f'event {i}: radial interval detached')
        if i>0 and c.predecessor_token!=cells[i-1].source_token:f.append(f'event {i}: literal predecessor detached')
    return list(dict.fromkeys(f))

def _identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def _zeros(n):return [[I(0) for _ in range(n)] for _ in range(n)]

def _smoke_image():
    st=JOINT._smoke_state()
    sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=JOINT.advance(st,sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='lineage')
    if not images:raise RuntimeError('joint estimator smoke emitted no image')
    return images[0]

def _smoke_lineage(radial:RadialCell):
    im=_smoke_image();n=18;P=_identity(n);x=[I(0) for _ in range(n)]
    e0=im.source_token+':e0';e1=im.source_token+':e1'
    c0=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=0,kind='prediction',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(0),radial_scale=radial.interval,event_source_token=e0,event_predecessor_token=im.predecessor_token)
    c1=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=1,kind='aw_floor',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.005),radial_scale=radial.interval,event_source_token=e1,event_predecessor_token=e0)
    return im,EventLineage(radial,im.source_token,(c0,c1))

def build():
    signed=SIGNED_MASTER.build();sf=SIGNED_MASTER.validate(signed)
    master=MASTER.build();mf=MASTER.validate(master)
    bias=BIASISS.build();bf=BIASISS.validate(bias)
    fp=FPISS.build();ff=FPISS.validate(fp)
    prefix=PREFIX.build();pf=PREFIX.validate(prefix)
    bad={k:v for k,v in (('signed_master',sf),('joint_master',mf),('BIAS1_ISS',bf),('binary32_ISS',ff),('prefix_bridge',pf)) if v}
    if bad:raise RuntimeError('lineage prerequisites failed: '+repr(bad))

    root=root_radial();p1=split_radial(root);p2=refine_partition(p1,p1[1].token)
    im,lineage=_smoke_lineage(root);lf=validate_event_lineage(lineage)
    if lf:raise RuntimeError('event-lineage smoke failed: '+repr(lf))
    F=_identity(18);Q=_zeros(18)
    cov=COV.prediction_transition_closed(lineage.cells[0],lineage.cells[1],F,Q)

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_signal_estimator_ancestry_required':True,
      'literal_event_ancestry_required':True,
      'continuous_radial_coordinate_retained':True,
      'radial_root_is_closed_unit_interval':root.interval.lo==0 and root.interval.hi==1,
      'binary_refinement_children_union_exact_parent':p1[0].interval.lo==0 and p1[0].interval.hi==p1[1].interval.lo and p1[1].interval.hi==1,
      'recursive_partition_exactly_covers_unit_interval':partition_covers_unit(p2),
      'radial_point_sampling_forbidden':True,
      'estimator_owned_event_lineage_smoke_pass':not lf,
      'reachable_covariance_event_continuity_smoke_pass':cov,
      'same_signal_period_sigma_relation_consumed':bool(im.state.statistics.wpe_velocity_C.lo>=0 and im.state.statistics.tuner_band_C.lo>=0),
      'real_exact_chord_signed_master_prerequisite_consumed':bool(signed['exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master']),
      'real_joint_sector_master_builder_consumed':bool(master['same_history_quadratic_graph_sector_assembler_available']),
      'physical_BIAS1_ISS_supply_consumed':bool(bias['same_w_enters_error_and_true_bias']),
      'binary32_Kalman_reset_ISS_consumed':bool(fp['additive_ISS_channel_complete_for_conditional_P4']),
      'synthetic_identity_master_may_promote':False,
      'production_same_history_event_master_maps_materialized_here':False,
      'all_admitted_BRMM_predecessor_cells_covered_here':False,
      'all_continuous_radial_cells_certified_here':False,
      'production_endpoint_augmented_LDLT_closed_here':False,
      'production_every_prefix_augmented_LDLT_closed_here':False,
      'production_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'materialize real D/M/J/B and chord/reset/projection sector maps for each estimator-owned literal event prefix over every BRMM predecessor and continuous radial partition cell; attach BIAS1 source_map and binary32 fp_map, then invoke the BRMM prefix LDLT bridge on every lineage'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_signal_estimator_ancestry_required','literal_event_ancestry_required','continuous_radial_coordinate_retained','radial_root_is_closed_unit_interval','binary_refinement_children_union_exact_parent','recursive_partition_exactly_covers_unit_interval','radial_point_sampling_forbidden','estimator_owned_event_lineage_smoke_pass','reachable_covariance_event_continuity_smoke_pass','same_signal_period_sigma_relation_consumed','real_exact_chord_signed_master_prerequisite_consumed','real_joint_sector_master_builder_consumed','physical_BIAS1_ISS_supply_consumed','binary32_Kalman_reset_ISS_consumed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('synthetic_identity_master_may_promote','production_same_history_event_master_maps_materialized_here','all_admitted_BRMM_predecessor_cells_covered_here','all_continuous_radial_cells_certified_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','production_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'radial_cover':d['recursive_partition_exactly_covers_unit_interval'],'event_lineage':d['estimator_owned_event_lineage_smoke_pass'],'real_master_maps':d['production_same_history_event_master_maps_materialized_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())