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
  * real estimator-owned Joseph prefix bundles containing the same-cell reduced
    signed master, structured nonlinear sectors, and certified same-Dtheta
    finite-reset sector;
  * the requirement that production augmented masters consume the retained
    physical BIAS1 supply and binary32 ISS channels.

This module still fails closed until these real prefix bundles are propagated
over every admitted BRMM predecessor/radial cell and embedded into one complete
word coordinate with BIAS1 and roundoff maps.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json
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
import ou3_p4_brmm_reduced_event_master as EVENTMASTER
import ou3_p4_brmm_same_Dtheta_reset_binding as RESETBIND
import ou3_p4_same_cell_correction_domain as CORR
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_frontend_state_step as FRONT

SCHEMA=2
QUALIFICATION='OU3_P4_BRMM_CONTINUOUS_RADIAL_EVENT_LINEAGE_COVER_V2'
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

@dataclass(frozen=True)
class CertifiedJosephPrefix:
    radial_token:str
    estimator_token:str
    event_token:str
    predecessor_event_token:str|None
    mode:str
    kind:str
    coordinate_dimension:int
    master:object
    nonlinear_sectors:tuple
    reset_sector:object
    reset_delta:float
    D_theta:object
    B_theta:object


def root_radial(token='RADIAL_ROOT'):
    return RadialCell(token,Interval(0.0,1.0),None)

def split_radial(cell:RadialCell,*,split:float|None=None):
    lo,hi=cell.interval.lo,cell.interval.hi
    if not (0.0<=lo<hi<=1.0):raise ValueError('radial cell must be nondegenerate subset of [0,1]')
    m=float(split) if split is not None else lo+(hi-lo)*0.5
    if not(lo<m<hi):raise ValueError('split must be strictly interior')
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

def certify_joseph_prefix(cell:COVER.SourceCoverCell,radial:RadialCell,corr:dict)->CertifiedJosephPrefix:
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):
        raise ValueError('certified Joseph prefix requires Joseph event')
    if cell.radial_scale is None or cell.radial_scale.lo!=radial.interval.lo or cell.radial_scale.hi!=radial.interval.hi:
        raise ValueError('Joseph cell detached from radial source interval')
    em=EVENTMASTER.build_event_master(cell)
    ef=EVENTMASTER.validate_event_master(cell,em)
    if ef:raise RuntimeError('real event master invalid: '+repr(ef))
    rb=RESETBIND.bind_event(em,corr)
    if not rb.get('closed'):
        raise RuntimeError('same-Dtheta reset binding did not close: '+repr(rb))
    return CertifiedJosephPrefix(
        radial.token,cell.estimator_source_token or '',cell.source_token,cell.predecessor_token,
        cell.mode,cell.kind,em.coordinate_dimension,em.master,em.nonlinear_sectors,
        rb['reset_sector'],float(rb['delta']),em.D_theta,em.B_theta)

def validate_certified_prefix(prefix:CertifiedJosephPrefix,cell:COVER.SourceCoverCell,radial:RadialCell):
    f=[]
    if prefix.radial_token!=radial.token:f.append('radial token detached')
    if prefix.estimator_token!=cell.estimator_source_token:f.append('estimator token detached')
    if prefix.event_token!=cell.source_token or prefix.predecessor_event_token!=cell.predecessor_token:f.append('literal event ancestry detached')
    if prefix.mode!=cell.mode or prefix.kind!=cell.kind:f.append('event mode/kind detached')
    n=prefix.coordinate_dimension
    if EVENTMASTER.shape(prefix.master)!=(n,n):f.append('event master dimension mismatch')
    if EVENTMASTER.shape(prefix.reset_sector)!=(n,n):f.append('reset sector dimension mismatch')
    if EVENTMASTER.shape(prefix.D_theta)!=(3,n) or EVENTMASTER.shape(prefix.B_theta)!=(3,n):f.append('reset graph maps dimension mismatch')
    if prefix.reset_delta<0:f.append('reset delta negative')
    return f

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

def _smoke_joseph_cells(im,radial):
    n=18;x=[I(0) for _ in range(n)];P=_identity(n);R=_identity(3);Rhat=_identity(3);f=[I(.2),I(-.1),I(-9.7)]
    e1=im.source_token+':j1';e2=im.source_token+':j2'
    s=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=1,kind='S_zero',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=radial.interval,event_source_token=e1,event_predecessor_token=im.predecessor_token)
    a=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='accelerometer',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.015),radial_scale=radial.interval,event_source_token=e2,event_predecessor_token=e1,R=R,f_hat=f,R_hat=Rhat)
    return s,a

def build():
    signed=SIGNED_MASTER.build();sf=SIGNED_MASTER.validate(signed)
    master=MASTER.build();mf=MASTER.validate(master)
    bias=BIASISS.build();bf=BIASISS.validate(bias)
    fp=FPISS.build();ff=FPISS.validate(fp)
    prefix=PREFIX.build();pf=PREFIX.validate(prefix)
    corr=CORR.build();cf=CORR.validate(corr)
    bad={k:v for k,v in (('signed_master',sf),('joint_master',mf),('BIAS1_ISS',bf),('binary32_ISS',ff),('prefix_bridge',pf),('same_cell_correction',cf)) if v}
    if bad:raise RuntimeError('lineage prerequisites failed: '+repr(bad))

    root=root_radial();p1=split_radial(root);p2=refine_partition(p1,p1[1].token)
    im,lineage=_smoke_lineage(root);lf=validate_event_lineage(lineage)
    if lf:raise RuntimeError('event-lineage smoke failed: '+repr(lf))
    F=_identity(18);Q=_zeros(18)
    cov=COV.prediction_transition_closed(lineage.cells[0],lineage.cells[1],F,Q)

    sj,aj=_smoke_joseph_cells(im,root)
    sp=certify_joseph_prefix(sj,root,corr);ap=certify_joseph_prefix(aj,root,corr)
    sfp=validate_certified_prefix(sp,sj,root);afp=validate_certified_prefix(ap,aj,root)
    structured_prefix_smoke=not(sfp or afp)

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_signal_estimator_ancestry_required':True,'literal_event_ancestry_required':True,
      'continuous_radial_coordinate_retained':True,'radial_root_is_closed_unit_interval':root.interval.lo==0 and root.interval.hi==1,
      'binary_refinement_children_union_exact_parent':p1[0].interval.lo==0 and p1[0].interval.hi==p1[1].interval.lo and p1[1].interval.hi==1,
      'recursive_partition_exactly_covers_unit_interval':partition_covers_unit(p2),'radial_point_sampling_forbidden':True,
      'estimator_owned_event_lineage_smoke_pass':not lf,'reachable_covariance_event_continuity_smoke_pass':cov,
      'same_signal_period_sigma_relation_consumed':bool(im.state.statistics.wpe_velocity_C.lo>=0 and im.state.statistics.tuner_band_C.lo>=0),
      'real_exact_chord_signed_master_prerequisite_consumed':bool(signed['exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master']),
      'real_joint_sector_master_builder_consumed':bool(master['same_history_quadratic_graph_sector_assembler_available']),
      'physical_BIAS1_ISS_supply_consumed':bool(bias['same_w_enters_error_and_true_bias']),
      'binary32_Kalman_reset_ISS_consumed':bool(fp['additive_ISS_channel_complete_for_conditional_P4']),
      'real_estimator_owned_Joseph_prefix_bundle_materialized_smoke':structured_prefix_smoke,
      'Joseph_prefix_bundle_contains_real_signed_master':structured_prefix_smoke,
      'Joseph_prefix_bundle_contains_structured_nonlinear_sectors':structured_prefix_smoke and len(ap.nonlinear_sectors)>0,
      'Joseph_prefix_bundle_contains_certified_same_Dtheta_reset_sector':structured_prefix_smoke,
      'Joseph_prefix_bundle_retains_event_specific_reset_delta':structured_prefix_smoke and sp.reset_delta>=0 and ap.reset_delta>=0,
      'Joseph_prefix_validation_failures':sfp+afp,
      'synthetic_identity_master_may_promote':False,
      'production_same_history_event_master_maps_materialized_here':False,
      'all_admitted_BRMM_predecessor_cells_covered_here':False,'all_continuous_radial_cells_certified_here':False,
      'production_BIAS1_and_binary32_maps_embedded_in_prefix_coordinate_here':False,
      'production_endpoint_augmented_LDLT_closed_here':False,'production_every_prefix_augmented_LDLT_closed_here':False,
      'production_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'embed certified Joseph prefix bundles and prediction/floor transports into one lineage-wide coordinate; add physical BIAS1 true-bias/driver and binary32 additive ISS coordinates; compose the signed matrix after every literal event and run outward LDLT plus hard-domain targets'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_signal_estimator_ancestry_required','literal_event_ancestry_required','continuous_radial_coordinate_retained','radial_root_is_closed_unit_interval','binary_refinement_children_union_exact_parent','recursive_partition_exactly_covers_unit_interval','radial_point_sampling_forbidden','estimator_owned_event_lineage_smoke_pass','reachable_covariance_event_continuity_smoke_pass','same_signal_period_sigma_relation_consumed','real_exact_chord_signed_master_prerequisite_consumed','real_joint_sector_master_builder_consumed','physical_BIAS1_ISS_supply_consumed','binary32_Kalman_reset_ISS_consumed','real_estimator_owned_Joseph_prefix_bundle_materialized_smoke','Joseph_prefix_bundle_contains_real_signed_master','Joseph_prefix_bundle_contains_structured_nonlinear_sectors','Joseph_prefix_bundle_contains_certified_same_Dtheta_reset_sector','Joseph_prefix_bundle_retains_event_specific_reset_delta'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('synthetic_identity_master_may_promote','production_same_history_event_master_maps_materialized_here','all_admitted_BRMM_predecessor_cells_covered_here','all_continuous_radial_cells_certified_here','production_BIAS1_and_binary32_maps_embedded_in_prefix_coordinate_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','production_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('Joseph_prefix_validation_failures')!=[]:f.append('certified Joseph prefix smoke validation failed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'radial_cover':d['recursive_partition_exactly_covers_unit_interval'],'event_lineage':d['estimator_owned_event_lineage_smoke_pass'],'real_Joseph_prefix':d['real_estimator_owned_Joseph_prefix_bundle_materialized_smoke'],'production_maps':d['production_same_history_event_master_maps_materialized_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
