#!/usr/bin/env python3
"""Same-graph S=0 attitude-correction domain certificate for P4.

For one estimator-owned S=0 Joseph event, form the actual same-cell correction
map

    d_theta = E_theta K H_S e = D_theta e

with K derived from that event's reachable P/H/R cell.  Prove

    ||d_theta|| <= delta h

from the declared hard-entry group balls by an affine S-procedure and outward
interval LDLT.  No independent/rowwise K bound, residual norm product, replay
point, or covariance-membership entry assumption is used.

This module closes the *per-cell* correction-domain implication for linear S=0
events.  It does not yet prove source-uniform coverage over every reachable
estimator/covariance predecessor cell, nor the accelerometer/vector chord case.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_mul,matrix_transpose
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_reset_signed_absorption_diagnostic as RESET
import ou3_p4_brmm_reduced_event_master as EVENT

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_SZERO_SAME_GRAPH_CORRECTION_DOMAIN_V1'

def I(x):return Interval.point(float(x))
def up(x):return math.nextafter(float(x),math.inf)

def _shape(A):return len(A),len(A[0]) if A else 0

def _subcolumns(A,cols):return [[row[j] for j in cols] for row in A]

def _gershgorin_upper(A):
    n,m=_shape(A)
    if n==0 or n!=m:raise ValueError('square matrix required')
    best=0.0
    for i in range(n):
        radius=sum(A[i][j].abs_upper() for j in range(n) if j!=i)
        best=max(best,up(A[i][i].hi+radius))
    return up(best)

def certify_szero_event(em:EVENT.EventMaster,margin:float=1.02):
    if em.kind!='S_zero':raise ValueError('S=0 event master required')
    n=18 if em.mode=='H' else 21
    if em.coordinate_layout[0] != ('e',0,n):raise ValueError('unexpected event state layout')
    if not (math.isfinite(margin) and margin>1.0):raise ValueError('strict margin >1 required')

    # Reduced implication coordinate zr=[e;h]. Reset-defect coordinates are not
    # part of the correction map and must not create artificial zero LDLT pivots.
    cols=list(range(n))+[em.h_index]
    Dr=_subcolumns(em.D_theta,cols)
    nr=n+1;hidx=n
    target=HARD.mapped_ball_target(Dr,hidx,1.0)  # rescaled below after delta chosen
    hard=HARD.hard_entry_iqcs('H18' if em.mode=='H' else 'A21',nr,h_index=hidx,state_offset=0)
    entry=ENTRY.build();bad=ENTRY.validate(entry)
    if bad:raise RuntimeError('hard-entry prerequisite failed: '+repr(bad))

    gram=matrix_mul(matrix_transpose(Dr),Dr)
    state_gram=[row[:n] for row in gram[:n]]
    L=_gershgorin_upper(state_gram)
    lam=up(max(L,1.0e-18)*margin)
    names=list(hard)
    sum_r2=0.0
    for name in names:
        r=float(entry['coordinate_radii'][name]);sum_r2=up(sum_r2+up(r*r))
    delta=up(math.sqrt(up(lam*sum_r2))*margin)

    # Rebuild the true target at the certified delta and use the exact declared
    # hard-ball IQCs. Equal multipliers are conservative but source-graph valid.
    target=HARD.mapped_ball_target(Dr,hidx,delta)
    multipliers=[lam for _ in names]
    ok,pivots=HARD.certify_strict_target(target,[hard[k] for k in names],multipliers)
    reset_domain_ok=False;reset_gain=None
    if ok:
        try:
            s=RESET.homogeneous_sector(float(entry['coordinate_radii']['attitude_cayley_norm']),delta)
            reset_domain_ok=bool(s['chart_safe'] and s['cayley_composition_denominator_lower']>0.0)
            reset_gain=float(s['reset_defect_over_correction_norm_upper'])
        except (ValueError,RuntimeError):
            reset_domain_ok=False
    return {
      'mode':em.mode,'event_source_token':em.source_token,'estimator_source_token':em.estimator_source_token,
      'same_cell_Dtheta_consumed':True,'rowwise_K_bound_used':False,'independent_residual_norm_used':False,
      'reduced_implication_coordinate':'[e;h]','hard_entry_groups':names,
      'Dtheta_gram_gershgorin_upper':L,'common_Sprocedure_multiplier':lam,
      'certified_delta':delta,'strict_outward_LDLT_closed':bool(ok),'LDLT_pivot_lowers':pivots,
      'parameterized_reset_sector_valid_at_certified_delta':reset_domain_ok,
      'reset_defect_over_correction_norm_upper':reset_gain,
    }

def build():
    im=EVENT._smoke_image();n=18;x=[I(0) for _ in range(n)];P=EVENT._identity(n)
    import ou3_p4_complete_brmm_source_cover_contract as COVER
    cell=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    em=EVENT.build_event_master(cell);ef=EVENT.validate_event_master(cell,em)
    if ef:raise RuntimeError('event master invalid: '+repr(ef))
    cert=certify_szero_event(em)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_P_H_R_K_event_consumed':True,'same_cell_Dtheta_target_certifier_materialized':True,
      'full_declared_hard_entry_balls_consumed':True,'covariance_membership_used':False,'rowwise_K_reset_domain_used':False,
      'smoke_certificate':cert,'smoke_same_graph_correction_domain_closed':cert['strict_outward_LDLT_closed'],
      'smoke_reset_sector_valid_at_certified_delta':cert['parameterized_reset_sector_valid_at_certified_delta'],
      'production_all_reachable_Szero_cells_certified_here':False,'production_accelerometer_vector_correction_domain_closed_here':False,
      'production_source_uniform_correction_domain_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'apply this same-Dtheta certifier to every reachable S=0 source cell; for accelerometer/vector retain chord/cross IQCs in the implication, then embed certified reset sectors into every literal prefix master'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_cell_P_H_R_K_event_consumed','same_cell_Dtheta_target_certifier_materialized','full_declared_hard_entry_balls_consumed','smoke_same_graph_correction_domain_closed','smoke_reset_sector_valid_at_certified_delta'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_membership_used','rowwise_K_reset_domain_used','production_all_reachable_Szero_cells_certified_here','production_accelerometer_vector_correction_domain_closed_here','production_source_uniform_correction_domain_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    c=d.get('smoke_certificate',{})
    if c.get('same_cell_Dtheta_consumed') is not True or c.get('rowwise_K_bound_used') is not False:f.append('smoke provenance invalid')
    if not(math.isfinite(float(c.get('certified_delta',math.nan))) and float(c['certified_delta'])>=0):f.append('certified delta invalid')
    if not c.get('LDLT_pivot_lowers') or min(map(float,c['LDLT_pivot_lowers']))<=0:f.append('strict LDLT pivot lost')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'same_graph':d['smoke_same_graph_correction_domain_closed'],'reset_valid':d['smoke_reset_sector_valid_at_certified_delta'],'delta':d['smoke_certificate']['certified_delta'],'source_uniform':d['production_source_uniform_correction_domain_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
