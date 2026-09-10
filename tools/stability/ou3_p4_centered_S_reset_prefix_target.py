#!/usr/bin/env python3
"""Centered S=0 same-cell reset-chart target for the production prefix proof.

Use augmented coordinate z=[h;e;truth/source;aux].  The pseudo residual selector
is Hc z=e_S-S_true=-S_hat.  For the same Joseph cell, d_theta=K_theta Hc z.
Rather than pre-bounding e_S or K independently, require the prefix S-procedure
to prove directly

    delta_*^2 h^2 - ||K_theta Hc z||^2 >= 0.

Here delta_* is chosen as the largest validated correction radius (to numerical
bisection tolerance) for which the deployed exact-reset homogeneous sector is
chart-safe over the full declared 45-degree Cayley state radius.  This radius
is a target/domain parameter, not an asserted source correction bound.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval,matrix_mul
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_reset_signed_absorption_diagnostic as RESETSECTOR
import ou3_p4_S_origin_literal_prefix_invariant as ORIGIN

SCHEMA=1
QUALIFICATION='OU3_P4_CENTERED_S_RESET_PREFIX_TARGET_V1'
S0=12

def _shape(A):return len(A),len(A[0]) if A else 0
def centered_selector(total_dim:int,error_offset:int,truth_offset:int,state_dim:int):
    if state_dim not in (18,21):raise ValueError('state dimension must be H18/A21')
    if total_dim<=max(error_offset+state_dim-1,truth_offset+state_dim-1):raise ValueError('augmented coordinate too small')
    z=Interval.point(0.0);one=Interval.point(1.0);minus=Interval.point(-1.0)
    H=[[z for _ in range(total_dim)] for _ in range(3)]
    for j in range(3):H[j][error_offset+S0+j]=one;H[j][truth_offset+S0+j]=minus
    return H
def correction_map(K,total_dim:int,error_offset:int,truth_offset:int,state_dim:int):
    r,c=_shape(K)
    if r!=state_dim or c!=3:raise ValueError('K must be state_dim x 3 from the same Joseph cell')
    H=centered_selector(total_dim,error_offset,truth_offset,state_dim)
    return matrix_mul([list(K[i]) for i in range(3)],H)
def chart_target(K,total_dim:int,h_index:int,error_offset:int,truth_offset:int,state_dim:int,delta_target:float):
    return HARD.mapped_ball_target(correction_map(K,total_dim,error_offset,truth_offset,state_dim),h_index,delta_target)
def _safe(q,d):
    try:
        s=RESETSECTOR.homogeneous_sector(q,d)
        return bool(s['chart_safe'] and s['homogeneous_sector_dominates_endpoint_absolute_bound'] and float(s['cayley_composition_denominator_lower'])>0)
    except (ValueError,RuntimeError,OverflowError):return False
def largest_safe_delta(q:float):
    lo=0.0;hi=3.0
    if not _safe(q,0.5):raise RuntimeError('known reset-sector smoke radius no longer safe')
    for _ in range(60):
        mid=(lo+hi)/2
        if _safe(q,mid):lo=mid
        else:hi=mid
    # Stay one representable float below the numerical boundary.
    d=math.nextafter(lo,0.0);sector=RESETSECTOR.homogeneous_sector(q,d)
    if not _safe(q,d):raise RuntimeError('selected reset target is not certified safe')
    return d,sector

def build():
    e=ENTRY.build();o=ORIGIN.build();bad={'entry':ENTRY.validate(e),'origin':ORIGIN.validate(o)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('centered S target prerequisites failed: '+repr(bad))
    q=float(e['coordinate_radii']['attitude_cayley_norm']);delta,sector=largest_safe_delta(q)
    # Point smoke on [h;e21;truth21].  The K value is arbitrary test data only;
    # production callers must pass K from their same P/H/R Joseph cell.
    n=43;I=Interval.point;K=[[I(0.0) for _ in range(3)] for _ in range(21)]
    for i in range(3):K[i][i]=I(0.1)
    H=centered_selector(n,1,22,21);D=correction_map(K,n,1,22,21);T=chart_target(K,n,0,1,22,21,delta)
    # Exact common-origin injection has equal S columns in e and truth; Hc kills it.
    origin_cancel=True
    for j in range(3):
        val=H[j][1+S0+j]+H[j][22+S0+j]
        origin_cancel=origin_cancel and val.lo==0.0 and val.hi==0.0
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'augmented_residual_selector':'Hc=[E_S,-E_S] on [e;truth/source]','same_cell_correction_map':'Dtheta=Etheta*K*Hc',
      'shared_origin_exactly_annihilated':origin_cancel,'literal_prefix_origin_invariant_consumed':o['finite_literal_prefix_composition_closed_by_induction'],
      'legacy_eS_radius_used_as_pseudo_residual':False,'independent_K_box_used':False,'scalar_S_residual_bound_required':False,
      'full_declared_attitude_cayley_norm':q,'validated_reset_chart_target_delta':delta,'validated_reset_sector_at_target':sector,
      'target_is_domain_goal_not_assumed_correction_bound':True,'same_graph_mapped_ball_target_constructor_available':len(T)==n and len(T[0])==n,
      'centered_selector_constructor_available':len(H)==3 and len(H[0])==n,'same_cell_Dtheta_constructor_available':len(D)==3 and len(D[0])==n,
      'production_S_zero_reset_chart_target_closed_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_PASS':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('shared_origin_exactly_annihilated','literal_prefix_origin_invariant_consumed','target_is_domain_goal_not_assumed_correction_bound','same_graph_mapped_ball_target_constructor_available','centered_selector_constructor_available','same_cell_Dtheta_constructor_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('legacy_eS_radius_used_as_pseudo_residual','independent_K_box_used','scalar_S_residual_bound_required','production_S_zero_reset_chart_target_closed_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    dlt=float(d.get('validated_reset_chart_target_delta',0));s=d.get('validated_reset_sector_at_target',{})
    if not(0.5<dlt<3.0):f.append('reset chart target delta invalid')
    if s.get('chart_safe') is not True or s.get('homogeneous_sector_dominates_endpoint_absolute_bound') is not True:f.append('reset sector target unsafe')
    return f
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'delta_target':d['validated_reset_chart_target_delta'],'mu':d['validated_reset_sector_at_target']['reset_defect_over_correction_norm_upper'],'origin_cancel':d['shared_origin_exactly_annihilated'],'production_target_closed':d['production_S_zero_reset_chart_target_closed_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
