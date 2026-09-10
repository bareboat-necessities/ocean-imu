#!/usr/bin/env python3
"""Bridge trusted complete-window kernel captures to the P4 event lineage.

The typed complete-window kernel is the shipping-order authority for Riccati
operations. With passive capture enabled it emits, from the SAME transition
call, literal event cells containing P_before/P_after and either F/Q or H/R.
This module binds those snapshots to the complete same-sample joint-estimator
image family and checks the literal shipping order

  prediction -> optional aw_floor -> due S_zero -> accelerometer -> async mag

including actual applied SpectralMSE R_S provenance.

The basic kernel frontend and the joint-estimator proof state are intentionally
not the same state machine: the latter carries central-moment/stillness proof
coordinates that the former does not. Therefore equality or whole-object
containment between their successor frontend objects is not a valid theorem
obligation. The actual common interface is the current active Riccati schedule
and the measurement geometry. Every retained joint image must bind to the same
captured current R_S, and the kernel must verify that it consumed that active
schedule before the current measurement.

This bridge intentionally does NOT fabricate a physical-error state from P.
Production SourceCoverCell construction remains fail-closed until the qualified
hard-entry physical-error set is propagated through the literal event sequence.
"""
from __future__ import annotations
import argparse,copy,json
from pathlib import Path

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_full_normal_live_word as WORD

SCHEMA=4
QUALIFICATION='OU3_P4_BRMM_TRUSTED_KERNEL_EVENT_CAPTURE_BRIDGE_V4'

def I(x):return Interval.point(float(x))
def _diag(n,v):
    z=I(0);d=I(v);return [[d if i==j else z for j in range(n)] for i in range(n)]
def _same_interval(a,b):return isinstance(a,Interval) and isinstance(b,Interval) and a.lo==b.lo and a.hi==b.hi
def _same_matrix(A,B):
    return len(A)==len(B) and all(len(A[i])==len(B[i]) and all(_same_interval(A[i][j],B[i][j]) for j in range(len(A[i]))) for i in range(len(A)))
def _sample():
    z=I(0)
    return KERNEL.SampleCoordinates(
      gyro_measurement=MAHONY.Vec3(I(.01),I(-.02),I(.005)),
      omega_body_corrected=(I(.01),I(-.02),I(.005)),
      specific_force=MAHONY.Vec3(I(.2),I(-.1),I(-9.75)),
      f_cog_body=(I(0),I(0),I(-9.80665)),
      R_wb=[[I(1),z,z],[z,I(1),z],[z,z,I(1)]],due_S=True,aw_floor_requested=True,
      magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20),I(0),I(40))),))

def capture_smoke():
    js=JOINT._smoke_state();sample=_sample();kc=KERNEL._process_constants()
    root=KERNEL.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=WORD.initialize_word('H',_diag(18,2)),A=WORD.initialize_word('A',_diag(21,2)),source_cell_id='root')
    succ,meta=KERNEL.advance_branch(root,sample,constants=kc,next_cell_prefix='capture',capture_riccati_event_cells=True)
    ji=JOINT.advance(js,FRONT.Sample(sample.gyro_measurement,sample.specific_force),gravity_ms2=kc.gravity,two_kp=kc.two_kp,two_ki=kc.two_ki,child_prefix='capture-joint')
    if not ji:raise RuntimeError('joint estimator emitted no same-sample image')
    if not succ:raise RuntimeError('trusted kernel emitted no same-sample successor')
    return sample,ji[0],meta,succ,ji

def schedule_binding(meta,images):
    h_s=next(c for c in meta['H_event_cells'] if c.kind=='S_zero')
    a_s=next(c for c in meta['A_event_cells'] if c.kind=='S_zero')
    per_image=[]
    for im in images:
        expected_h=WORD.R_S_zero(im.actual_rs_std_xyz_for_current_riccati)
        expected_a=WORD.R_S_zero(im.actual_rs_std_xyz_for_current_riccati)
        per_image.append({
          'source_token':im.source_token,
          'H18_current_RS_matches_captured_kernel_event':_same_matrix(h_s.R,expected_h),
          'A21_current_RS_matches_captured_kernel_event':_same_matrix(a_s.R,expected_a)})
    all_bound=bool(per_image) and all(x['H18_current_RS_matches_captured_kernel_event'] and x['A21_current_RS_matches_captured_kernel_event'] for x in per_image)
    return {
      'binding_relation':'joint current-active schedule -> exact captured kernel S=0 R; auxiliary successor frontend state is not equated',
      'kernel_verified_same_active_schedule_before_measurement':bool(meta['same_active_schedule_verified']),
      'kernel_verified_same_actual_RS_before_measurement':bool(meta['same_actual_RS_verified']),
      'every_joint_image_current_schedule_bound_to_kernel_Riccati_event':all_bound,
      'joint_image_bindings':per_image,
      'kernel_successor_family_nonempty':True,
      'joint_estimator_image_family_nonempty':True,
      'whole_frontend_successor_equality_required':False,
      'whole_frontend_successor_containment_required':False}
def summarize(cells,images,sample,mode):
    kinds=[c.kind for c in cells]
    expected=['prediction','aw_floor','S_zero','accelerometer','magnetometer']
    if kinds!=expected:raise RuntimeError(f'{mode} literal capture order changed: {kinds}')
    n=18 if mode=='H' else 21;chain=True
    for i,c in enumerate(cells):
        if len(c.P_before)!=n or len(c.P_after)!=n:chain=False
        if i and not _same_matrix(c.P_before,cells[i-1].P_after):chain=False
    s=next(c for c in cells if c.kind=='S_zero');acc=next(c for c in cells if c.kind=='accelerometer');mag=next(c for c in cells if c.kind=='magnetometer');pred=next(c for c in cells if c.kind=='prediction')
    rs_match=all(_same_matrix(s.R,WORD.R_S_zero(im.actual_rs_std_xyz_for_current_riccati)) for im in images)
    return {
      'mode':mode,'literal_kinds':kinds,'covariance_chain_contiguous':chain,
      'prediction_F_Q_captured':pred.F is not None and pred.Q is not None,
      'S_zero_H_R_captured':s.H is not None and s.R is not None,'S_zero_actual_RS_flag':s.actual_rs_from_committed_schedule,
      'S_zero_R_matches_every_joint_image_current_active_schedule':rs_match,
      'accelerometer_H_R_captured':acc.H is not None and acc.R is not None,
      'accelerometer_geometry_matches_sample':_same_matrix(acc.H,WORD.H_accelerometer(mode,sample.f_cog_body,sample.R_wb)),
      'magnetometer_H_R_captured':mag.H is not None and mag.R is not None,'magnetometer_index':mag.magnetic_event_index}
def build():
    sample,_,meta,succ,images=capture_smoke();binding=schedule_binding(meta,images);H=summarize(meta['H_event_cells'],images,sample,'H');A=summarize(meta['A_event_cells'],images,sample,'A')
    all_ok=all(x for d in (H,A) for k,x in d.items() if k not in ('mode','literal_kinds','magnetometer_index'))
    all_bound=binding['every_joint_image_current_schedule_bound_to_kernel_Riccati_event'] and binding['kernel_verified_same_active_schedule_before_measurement'] and binding['kernel_verified_same_actual_RS_before_measurement']
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_kernel_passive_capture_consumed':True,'same_transition_call_not_replay':True,'shipping_literal_event_order_captured':True,'H18_A21_same_sample_capture':True,
      'same_sample_joint_estimator_image_family_materialized':True,'joint_image_branch_count':len(images),'kernel_frontend_successor_count':len(succ),
      'same_sample_schedule_binding':binding,'all_joint_images_bound_to_current_kernel_Riccati_schedule':all_bound,
      'frontend_auxiliary_state_not_used_as_cross_machine_identity':True,
      'favorable_joint_estimator_successor_selected_for_production':False,
      'current_active_schedule_matches_kernel_Riccati_schedule':bool(meta['same_active_schedule_verified']),'current_actual_RS_matches_kernel_Riccati_schedule':bool(meta['same_actual_RS_verified']),
      'H18':H,'A21':A,'smoke_capture_all_same_cell_objects_closed':all_ok,
      'covariance_used_as_physical_error_membership':False,'hard_entry_state_transport_fabricated_here':False,
      'production_joint_estimator_all_successors_bound_here':False,'production_physical_error_event_states_propagated_here':False,
      'production_SourceCoverCells_emitted_here':False,'production_complete_event_lineage_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'propagate the qualified hard-entry physical-error enclosure through these exact captured operations for every retained joint image; construct estimator-owned SourceCoverCells with captured P/H/R and geometry, then close radial/source continuation rather than selecting a successor'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_passive_capture_consumed','same_transition_call_not_replay','shipping_literal_event_order_captured','H18_A21_same_sample_capture','same_sample_joint_estimator_image_family_materialized','all_joint_images_bound_to_current_kernel_Riccati_schedule','frontend_auxiliary_state_not_used_as_cross_machine_identity','current_active_schedule_matches_kernel_Riccati_schedule','current_actual_RS_matches_kernel_Riccati_schedule','smoke_capture_all_same_cell_objects_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    b=d.get('same_sample_schedule_binding',{})
    for k in ('kernel_verified_same_active_schedule_before_measurement','kernel_verified_same_actual_RS_before_measurement','every_joint_image_current_schedule_bound_to_kernel_Riccati_event','kernel_successor_family_nonempty','joint_estimator_image_family_nonempty'):
        if b.get(k) is not True:f.append('binding '+k+' not true')
    for k in ('whole_frontend_successor_equality_required','whole_frontend_successor_containment_required'):
        if b.get(k) is not False:f.append('binding '+k+' not false')
    for mode in ('H18','A21'):
        m=d.get(mode,{})
        for k in ('covariance_chain_contiguous','prediction_F_Q_captured','S_zero_H_R_captured','S_zero_actual_RS_flag','S_zero_R_matches_every_joint_image_current_active_schedule','accelerometer_H_R_captured','accelerometer_geometry_matches_sample','magnetometer_H_R_captured'):
            if m.get(k) is not True:f.append(mode+' '+k+' not true')
    for k in ('favorable_joint_estimator_successor_selected_for_production','covariance_used_as_physical_error_membership','hard_entry_state_transport_fabricated_here','production_joint_estimator_all_successors_bound_here','production_physical_error_event_states_propagated_here','production_SourceCoverCells_emitted_here','production_complete_event_lineage_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'capture':d['smoke_capture_all_same_cell_objects_closed'],'schedule_bound':d['all_joint_images_bound_to_current_kernel_Riccati_schedule'],'joint_images':d['joint_image_branch_count'],'kernel_successors':d['kernel_frontend_successor_count'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
