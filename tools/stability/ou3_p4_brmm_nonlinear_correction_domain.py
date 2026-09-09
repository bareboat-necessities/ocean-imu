#!/usr/bin/env python3
"""Same-graph accelerometer/vector correction-domain certificate for P4.

For an estimator-owned nonlinear Joseph event, retain the exact event coordinate
from ``ou3_p4_brmm_reduced_event_master``.  The attitude correction is

    d_theta = D_theta z,

with D_theta built from the SAME P/H/R/K cell and the structured residual map.
We seek a strict affine S-procedure certificate

    delta^2 h^2 - ||D_theta z||^2
      - sum_i lambda_i Pi_i  > 0,      lambda_i >= 0,

where the premise IQCs include the full declared hard-entry balls and *all*
nonlinear chord/cross-product graph sectors carried by the event master.

This is deliberately stronger and more faithful than multiplying a residual
norm by a rowwise K norm.  Equality graph constraints are retained in their
negative-Gram form; chord retention and cross-product orthogonality/norm sectors
remain on the same augmented coordinate.  Candidate multipliers are searched
only to discover a strict outward-LDLT certificate.  Failure to find one is
reported as an obstruction and never converted into a proof claim.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_brmm_reduced_event_master as EVENT
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_NONLINEAR_SAME_GRAPH_CORRECTION_DOMAIN_V1'

def _groups(em):
    hard=HARD.hard_entry_iqcs('H18' if em.mode=='H' else 'A21',em.coordinate_dimension,h_index=em.h_index,state_offset=0)
    nonlinear=dict(em.nonlinear_sectors)
    eq=[Pi for name,Pi in nonlinear.items() if name.startswith('p_equals_') or name.endswith('_plus') or name.endswith('_minus')]
    chord=[Pi for name,Pi in nonlinear.items() if name.startswith('chord_') and Pi not in eq]
    cross=[Pi for name,Pi in nonlinear.items() if name.startswith('cross_') and Pi not in eq]
    return hard,nonlinear,eq,chord,cross

def _candidate_multipliers(hard,nonlinear,lh,le,lc,lx):
    vals=[]
    for _ in hard:vals.append(lh)
    for name in nonlinear:
        if name.startswith('p_equals_') or name.endswith('_plus') or name.endswith('_minus'):vals.append(le)
        elif name.startswith('chord_'):vals.append(lc)
        elif name.startswith('cross_'):vals.append(lx)
        else:vals.append(lc)
    return vals

def certify_event(em,*,delta_candidates=None,multiplier_grid=None):
    if em.kind not in ('accelerometer','magnetometer'):raise ValueError('nonlinear Joseph event required')
    hard,nonlinear,eq,chord,cross=_groups(em)
    premises=list(hard.values())+list(nonlinear.values())
    if not nonlinear:raise ValueError('nonlinear event lost chord/cross sectors')
    if delta_candidates is None:
        # Stay entirely inside the exact-reset chart utility.  Values are not
        # assumed certificates; every candidate is checked by strict LDLT.
        cap=float(RESET.CAYLEY_MONOTONE_NORM_MAX)
        delta_candidates=[cap*x for x in (.20,.30,.40,.50,.60,.70,.80,.90,.98)]
    if multiplier_grid is None:
        multiplier_grid=(1e-4,1e-2,1.0,1e2,1e4)
    best=None;attempts=0
    for delta in delta_candidates:
        if not(0.0<float(delta)<=float(RESET.CAYLEY_MONOTONE_NORM_MAX)):continue
        target=HARD.mapped_ball_target(em.D_theta,em.h_index,float(delta))
        # Structured four-scale search.  Equal scale within a semantic IQC
        # family keeps the search deterministic and source-independent.
        for lh in multiplier_grid:
          for le in multiplier_grid:
            for lc in multiplier_grid:
              for lx in multiplier_grid:
                mult=_candidate_multipliers(hard,nonlinear,lh,le,lc,lx);attempts+=1
                ok,piv=HARD.certify_strict_target(target,premises,mult)
                if ok:
                    best={'delta':float(delta),'hard_multiplier':lh,'equality_multiplier':le,
                          'chord_multiplier':lc,'cross_multiplier':lx,'pivot_lowers':piv,
                          'minimum_pivot_lower':min(piv),'attempt':attempts}
                    return {'closed':True,'certificate':best,'attempts':attempts,
                            'hard_premise_count':len(hard),'nonlinear_premise_count':len(nonlinear),
                            'equality_like_count':len(eq),'chord_count':len(chord),'cross_count':len(cross)}
    return {'closed':False,'certificate':None,'attempts':attempts,
            'hard_premise_count':len(hard),'nonlinear_premise_count':len(nonlinear),
            'equality_like_count':len(eq),'chord_count':len(chord),'cross_count':len(cross)}

def _mag_smoke_cell(im):
    n=18;x=[EVENT.I(0) for _ in range(n)];P=EVENT._identity(n);R=EVENT._identity(3);m=[EVENT.I(.3),EVENT.I(-.1),EVENT.I(.2)]
    return EVENT.COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=4,kind='magnetometer',state=x,P=P,dt_s=EVENT.I(.005),pseudo_elapsed_s=EVENT.I(.02),radial_scale=EVENT.Interval(0,1),event_source_token=im.source_token+':e4',event_predecessor_token=im.source_token+':e3',R=R,m_body=m)
def build():
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    im=EVENT._smoke_image();acc=EVENT.build_event_master(EVENT._accel_smoke_cell(im));mag=EVENT.build_event_master(_mag_smoke_cell(im))
    af=EVENT.validate_event_master(EVENT._accel_smoke_cell(im),acc);mf=EVENT.validate_event_master(_mag_smoke_cell(im),mag)
    if af or mf:raise RuntimeError(f'event master invalid accel={af} mag={mf}')
    a=certify_event(acc);m=certify_event(mag)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_Dtheta_consumed':True,'full_declared_hard_entry_balls_consumed':True,
      'exact_chord_sectors_retained_in_correction_implication':True,
      'accelerometer_c_cross_aw_sectors_retained_in_correction_implication':True,
      'rowwise_K_bound_used':False,'independent_residual_norm_used':False,'covariance_membership_used':False,
      'strict_outward_LDLT_is_only_success_condition':True,'multiplier_search_is_diagnostic_not_assumption':True,
      'accelerometer':a,'magnetometer':m,
      'smoke_accelerometer_same_graph_closed':bool(a['closed']),
      'smoke_magnetometer_same_graph_closed':bool(m['closed']),
      'production_all_reachable_nonlinear_cells_certified_here':False,
      'production_source_uniform_nonlinear_correction_domain_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':('if smoke certificates close, apply the same certifier to every reachable estimator-owned nonlinear event cell and radial segment; '
        'if either smoke fails, refine the multiplier family or retain additional exact graph equalities before any reset-domain promotion')}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_cell_Dtheta_consumed','full_declared_hard_entry_balls_consumed','exact_chord_sectors_retained_in_correction_implication','accelerometer_c_cross_aw_sectors_retained_in_correction_implication','strict_outward_LDLT_is_only_success_condition','multiplier_search_is_diagnostic_not_assumption'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_bound_used','independent_residual_norm_used','covariance_membership_used','production_all_reachable_nonlinear_cells_certified_here','production_source_uniform_nonlinear_correction_domain_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    # Smoke closure is intentionally not a schema-validity requirement. A
    # genuine failure is evidence about the current S-procedure, not malformed proof code.
    for name in ('accelerometer','magnetometer'):
        r=d.get(name,{})
        if int(r.get('attempts',0))<=0:f.append(name+' search did not execute')
        if r.get('closed'):
            c=r.get('certificate') or {}
            if not c.get('pivot_lowers') or min(map(float,c['pivot_lowers']))<=0:f.append(name+' claimed closure without positive outward pivots')
            if not(0<float(c.get('delta',0))<=float(RESET.CAYLEY_MONOTONE_NORM_MAX)):f.append(name+' certified delta outside reset chart')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'accel_closed':d['smoke_accelerometer_same_graph_closed'],'mag_closed':d['smoke_magnetometer_same_graph_closed'],'accel_attempts':d['accelerometer']['attempts'],'mag_attempts':d['magnetometer']['attempts'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
