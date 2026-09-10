#!/usr/bin/env python3
"""A21 real Joseph prefix augmented with physical beta and active projection graph.

Start from one estimator-owned A21 Joseph event master. Extend its common
coordinate by beta_true(3) and projected f_ba(3). The pre-projection bias error
uses the SAME event correction d=Kq. The deployed radial projection graph, true
bias hard bound, nonlinear chord sectors, and a finite-reset sector all share
one coordinate.

The reset chart is now handled correctly as a first-exit obligation: this module
embeds the SAME-Dtheta target ||Dtheta z||<=3 h and the reset IQC that is valid
under that target. It does not pre-assert the target through the obsolete huge
scalar covariance/residual correction ceiling. The common production augmented
master must prove the target before it may consume the reset sector.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_sub
import ou3_p4_brmm_reduced_event_master as EVENT
import ou3_p4_brmm_same_Dtheta_reset_binding as RESETBIND
import ou3_p4_brmm_projection_graph_iqc as PROJIQC
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_bias1_family as BIAS1

SCHEMA=2
QUALIFICATION='OU3_P4_BRMM_A21_REAL_PROJECTION_PREFIX_V2'

def I(x):return Interval.point(float(x))
def shape(A):return len(A),len(A[0]) if A else 0
def zeros(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def embed_square(A,n):
    r,c=shape(A)
    if r!=c or r>n:raise ValueError('square embedding mismatch')
    out=zeros(n,n)
    for i in range(r):
        for j in range(c):out[i][j]=A[i][j]
    return out
def extend_map(A,n):
    r,c=shape(A)
    if c>n:raise ValueError('map embedding mismatch')
    return [list(row)+[I(0) for _ in range(n-c)] for row in A]
def selector3(n,offset):
    if offset<0 or offset+3>n:raise ValueError('selector outside coordinate')
    return [[I(1 if j==offset+i else 0) for j in range(n)] for i in range(3)]

def build_augmented_prefix(cell):
    if cell.mode!='A':raise ValueError('A21 source cell required')
    em=EVENT.build_event_master(cell);ef=EVENT.validate_event_master(cell,em)
    if ef:raise RuntimeError('base A21 event master invalid: '+repr(ef))
    rb=RESETBIND.bind_event_first_exit(em)
    if not rb.get('closed'):raise RuntimeError('same-Dtheta first-exit reset structure failed: '+repr(rb))
    n0=em.coordinate_dimension;beta0=n0;f0=n0+3;n=n0+6
    master=embed_square(em.master,n)
    nonlinear=tuple((name,embed_square(Pi,n)) for name,Pi in em.nonlinear_sectors)
    reset_sector=embed_square(rb['reset_sector'],n)
    reset_target=embed_square(rb['correction_domain_target'],n)
    D=extend_map(em.D,n);E=extend_map(em.E,n)
    Eba=[list(row) for row in E[18:21]];Dba=[list(row) for row in D[18:21]]
    preba=matrix_sub(Eba,Dba)
    Beta=selector3(n,beta0);F=selector3(n,f0)
    projection_sector=PROJIQC.projection_iqc(preba,Beta,F)
    fam=BIAS1.build();bf=BIAS1.validate(fam)
    if bf:raise RuntimeError('BIAS1 prerequisite failed: '+repr(bf))
    beta_bound=HARD.ball_iqc(n,em.h_index,(beta0,beta0+1,beta0+2),float(fam['true_bias_norm_upper_mps2']))
    return {'base_event_master':master,'nonlinear_sectors':nonlinear,'reset_sector':reset_sector,
            'reset_correction_domain_target':reset_target,'projection_sector':projection_sector,'beta_bound_sector':beta_bound,
            'coordinate_dimension':n,'base_dimension':n0,'beta_offset':beta0,'projected_ba_offset':f0,
            'h_index':em.h_index,'preprojection_ba_map':preba,'beta_map':Beta,'projected_ba_map':F,
            'reset_utility_delta':float(rb['delta']),'event_token':cell.source_token,
            'estimator_token':cell.estimator_source_token,'radial_scale':cell.radial_scale,
            'reset_correction_domain_target_proved_here':False}

def _smoke_cell():
    im=EVENT._smoke_image();n=21;x=[I(0) for _ in range(n)];P=EVENT._identity(n);R=EVENT._identity(3);Rhat=EVENT._identity(3);f=[I(.2),I(-.1),I(-9.7)];beta=[I(.02),I(-.01),I(.015)]
    return EVENT.COVER.source_cell_from_joint_image(im,mode='A',sample_index=0,event_ordinal=3,kind='accelerometer',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.015),radial_scale=Interval(0,1),event_source_token=im.source_token+':e3',event_predecessor_token=im.source_token+':e2',R=R,f_hat=f,R_hat=Rhat,true_bias=beta,bias_projection_limit=.4)
def build():
    proj=PROJIQC.build();pf=PROJIQC.validate(proj)
    if pf:raise RuntimeError('projection prefix prerequisite failed: '+repr(pf))
    cell=_smoke_cell();a=build_augmented_prefix(cell);n=a['coordinate_dimension']
    matrices=[a['base_event_master'],a['reset_sector'],a['reset_correction_domain_target'],a['projection_sector'],a['beta_bound_sector']]+[Pi for _,Pi in a['nonlinear_sectors']]
    dims=all(shape(M)==(n,n) for M in matrices)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'estimator_owned_A21_event_consumed':a['estimator_token']==cell.estimator_source_token,
      'same_event_D_equals_Kq_bias_correction_consumed':True,'preprojection_ba_is_e_ba_minus_same_cell_d_ba':True,
      'physical_BIAS1_beta_coordinate_materialized':True,'projected_ba_output_coordinate_materialized':True,
      'global_active_projection_IQC_attached':True,'saturated_unsaturated_Clarke_branches_all_retained':True,
      'BIAS1_true_bias_hard_bound_attached_on_same_h':True,
      'same_Dtheta_reset_sector_embedded':True,'same_Dtheta_first_exit_correction_target_embedded':True,
      'reset_target_preproved_by_scalar_ceiling':False,'reset_target_must_be_closed_by_common_augmented_master':True,
      'structured_chord_cross_sectors_preserved_under_embedding':len(a['nonlinear_sectors'])>0,
      'all_augmented_matrices_share_one_coordinate':dims,'augmented_coordinate_dimension':n,
      'reset_utility_delta':a['reset_utility_delta'],'inactive_projection_assumed':False,'rowwise_K_bound_used':False,
      'independent_beta_slots_per_event_used':False,'production_BIAS1_driver_recurrence_embedded_here':False,
      'production_binary32_forcing_embedded_here':False,'production_complete_prefix_storage_delta_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'carry beta_true and BIAS family driver state across prediction prefixes, add binary32 ISS forcing, and prove the embedded same-Dtheta correction target plus endpoint/every-prefix storage by the common augmented LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('estimator_owned_A21_event_consumed','same_event_D_equals_Kq_bias_correction_consumed','preprojection_ba_is_e_ba_minus_same_cell_d_ba','physical_BIAS1_beta_coordinate_materialized','projected_ba_output_coordinate_materialized','global_active_projection_IQC_attached','saturated_unsaturated_Clarke_branches_all_retained','BIAS1_true_bias_hard_bound_attached_on_same_h','same_Dtheta_reset_sector_embedded','same_Dtheta_first_exit_correction_target_embedded','reset_target_must_be_closed_by_common_augmented_master','structured_chord_cross_sectors_preserved_under_embedding','all_augmented_matrices_share_one_coordinate'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('reset_target_preproved_by_scalar_ceiling','inactive_projection_assumed','rowwise_K_bound_used','independent_beta_slots_per_event_used','production_BIAS1_driver_recurrence_embedded_here','production_binary32_forcing_embedded_here','production_complete_prefix_storage_delta_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not(math.isfinite(float(d.get('reset_utility_delta',math.nan))) and 0<float(d['reset_utility_delta'])<=3):f.append('reset utility delta invalid')
    if int(d.get('augmented_coordinate_dimension',0))<=0:f.append('coordinate dimension invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'A21_projection_prefix':d['all_augmented_matrices_share_one_coordinate'],'dim':d['augmented_coordinate_dimension'],'delta':d['reset_utility_delta'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
