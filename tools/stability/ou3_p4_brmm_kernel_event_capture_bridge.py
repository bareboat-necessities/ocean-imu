#!/usr/bin/env python3
"""Bridge trusted complete-window kernel captures to the P4 event lineage.

The typed complete-window kernel is the shipping-order authority for Riccati
operations.  With passive capture enabled it emits, from the SAME transition
call, literal event cells containing P_before/P_after and either F/Q or H/R.
This module binds those snapshots to the same-sample joint estimator Image and
checks that

  prediction -> optional aw_floor -> due S_zero -> accelerometer -> async mag

is preserved exactly, including actual applied SpectralMSE R_S provenance.

This bridge intentionally does NOT fabricate a physical-error state from P.
Production SourceCoverCell construction therefore remains fail-closed until the
qualified hard-entry physical-error set is propagated through the literal event
sequence.  Covariance consistency is never used as an entry-membership proxy.
"""
from __future__ import annotations
import argparse,copy,json
from pathlib import Path

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_shipping_prediction_primitives as PRED
import ou3_brmm_full_normal_live_word as WORD

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_TRUSTED_KERNEL_EVENT_CAPTURE_BRIDGE_V1'

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
    # Every joint image shares the current committed schedule. Choose one only
    # for interface smoke; no production branch is discarded or certified here.
    image=ji[0]
    return sample,image,meta,succ,ji

def summarize(cells,image,sample,mode):
    kinds=[c.kind for c in cells]
    expected=['prediction','aw_floor','S_zero','accelerometer','magnetometer']
    if kinds!=expected:raise RuntimeError(f'{mode} literal capture order changed: {kinds}')
    n=18 if mode=='H' else 21
    chain=True
    for i,c in enumerate(cells):
        if len(c.P_before)!=n or len(c.P_after)!=n:chain=False
        if i and not _same_matrix(c.P_before,cells[i-1].P_after):chain=False
    s=next(c for c in cells if c.kind=='S_zero')
    acc=next(c for c in cells if c.kind=='accelerometer')
    mag=next(c for c in cells if c.kind=='magnetometer')
    pred=next(c for c in cells if c.kind=='prediction')
    rs_expected=WORD.R_S_zero(image.actual_rs_std_xyz_for_current_riccati)
    return {
      'mode':mode,'literal_kinds':kinds,'covariance_chain_contiguous':chain,
      'prediction_F_Q_captured':pred.F is not None and pred.Q is not None,
      'S_zero_H_R_captured':s.H is not None and s.R is not None,
      'S_zero_actual_RS_flag':s.actual_rs_from_committed_schedule,
      'S_zero_R_matches_joint_image_active_schedule':_same_matrix(s.R,rs_expected),
      'accelerometer_H_R_captured':acc.H is not None and acc.R is not None,
      'accelerometer_geometry_matches_sample':_same_matrix(acc.H,WORD.H_accelerometer(mode,sample.f_cog_body,sample.R_wb)),
      'magnetometer_H_R_captured':mag.H is not None and mag.R is not None,
      'magnetometer_index':mag.magnetic_event_index,
    }
def build():
    sample,image,meta,succ,images=capture_smoke();H=summarize(meta['H_event_cells'],image,sample,'H');A=summarize(meta['A_event_cells'],image,sample,'A')
    all_ok=all(x for d in (H,A) for k,x in d.items() if k not in ('mode','literal_kinds','magnetometer_index'))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_kernel_passive_capture_consumed':True,'same_transition_call_not_replay':True,
      'shipping_literal_event_order_captured':True,'H18_A21_same_sample_capture':True,
      'same_sample_joint_estimator_image_materialized':True,'joint_image_branch_count':len(images),
      'favorable_joint_estimator_successor_selected_for_production':False,
      'current_active_schedule_matches_kernel_Riccati_schedule':bool(meta['same_active_schedule_verified']),
      'current_actual_RS_matches_kernel_Riccati_schedule':bool(meta['same_actual_RS_verified']),
      'H18':H,'A21':A,'smoke_capture_all_same_cell_objects_closed':all_ok,
      'covariance_used_as_physical_error_membership':False,'hard_entry_state_transport_fabricated_here':False,
      'production_joint_estimator_all_successors_bound_here':False,
      'production_physical_error_event_states_propagated_here':False,
      'production_SourceCoverCells_emitted_here':False,'production_complete_event_lineage_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'propagate the qualified hard-entry physical-error enclosure through these exact captured operations; for each Joseph boundary construct an estimator-owned SourceCoverCell with P_before plus captured H/R and source geometry, retaining every joint estimator successor rather than selecting the smoke image'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_passive_capture_consumed','same_transition_call_not_replay','shipping_literal_event_order_captured','H18_A21_same_sample_capture','same_sample_joint_estimator_image_materialized','current_active_schedule_matches_kernel_Riccati_schedule','current_actual_RS_matches_kernel_Riccati_schedule','smoke_capture_all_same_cell_objects_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for mode in ('H18','A21'):
        m=d.get(mode,{})
        for k in ('covariance_chain_contiguous','prediction_F_Q_captured','S_zero_H_R_captured','S_zero_actual_RS_flag','S_zero_R_matches_joint_image_active_schedule','accelerometer_H_R_captured','accelerometer_geometry_matches_sample','magnetometer_H_R_captured'):
            if m.get(k) is not True:f.append(mode+' '+k+' not true')
    for k in ('favorable_joint_estimator_successor_selected_for_production','covariance_used_as_physical_error_membership','hard_entry_state_transport_fabricated_here','production_joint_estimator_all_successors_bound_here','production_physical_error_event_states_propagated_here','production_SourceCoverCells_emitted_here','production_complete_event_lineage_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'capture':d['smoke_capture_all_same_cell_objects_closed'],'H':d['H18']['literal_kinds'],'A':d['A21']['literal_kinds'],'state_transport':d['production_physical_error_event_states_propagated_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
