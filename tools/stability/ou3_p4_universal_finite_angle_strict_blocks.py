#!/usr/bin/env python3
"""Materialize the universal finite-angle H18/A21 strict matrix blocks for P4.

The existing finite-angle producers prove strict full-matrix margins but retain
only summary pivots.  The nonlinear augmented S-procedure needs the actual
outward matrices as its top-left homogeneous state block.  This module rebuilds
*the same* H18 cell matrix used by ``ou3_brmm_h18_prior_free_completion`` with
the already-certified finite-angle completion penalty, then appends the exact
first-active A21 bias direct-sum margin.

No new source family, scalar contraction surrogate, trajectory replay, domain
shrink, or alternate storage is introduced.  The matrices are theorem-facing
API objects; the CLI writes only scalar metadata because Interval matrices are
not JSON artifacts.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_brmm_h18_prior_free_completion as HPF
import ou3_brmm_riccati_tube as TUBE
import ou3_brmm_riccati_tube_factored as FACTORED
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_full_process_ucc as PROCESS
import ou3_p4_finite_angle_h18_universal as H18FA
import ou3_p4_finite_angle_a21_universal as A21FA

QUALIFICATION='OU3_P4_UNIVERSAL_FINITE_ANGLE_STRICT_BLOCKS_V1'


def I(x:float)->Interval:
    return Interval.outward_bounds(float(x),float(x))


def _zero(n:int):
    z=I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _h18_matrix(x:Interval, *, process:dict, dynamic:dict, penalty_physical:float):
    """Exact matrix body of HPF._full_H18_cell, returned instead of discarded."""
    delta=HPF.USEFUL_GATE
    one_minus=I(math.nextafter(1.0-delta,-math.inf))
    M=_zero(18)
    q_att=float(process['attitude_gyro_bias']['Q_attitude_gyro_bias_lambda_min_lower'])
    att_diag=math.nextafter((1.0-delta)*q_att-float(penalty_physical),-math.inf)
    if not att_diag>0.0:raise RuntimeError('finite-angle strict block lost attitude margin')
    for i in range(6):M[i][i]=I(att_diag)
    inv=dynamic['dynamic_invariant']
    sigma_floor=float(inv['sigma_aw_filter_mps2'][0])
    h=float(dynamic['validated_rate_and_jump_bounds']['dt_s'])
    if not(sigma_floor>0.0 and h>0.0):raise RuntimeError('invalid BRMM process scale')
    scales=[sigma_floor*h,sigma_floor*h*h,sigma_floor*h*h*h,sigma_floor]
    B=FACTORED.step_scaled_q_over_x(x);qscale=I(x.lo)
    Maxis=[[one_minus*qscale*B[i][j] for j in range(4)] for i in range(4)]
    for i in range(4):
        denom=TUBE.down(scales[i]*scales[i])
        if not denom>0.0:raise RuntimeError('finite-angle congruence scale lost positivity')
        Maxis[i][i]=Maxis[i][i]-I(TUBE.up(float(penalty_physical)/denom))
    idx=(HPF.OFF_V,HPF.OFF_P,HPF.OFF_S,HPF.OFF_AW)
    for axis in range(3):
        for i in range(4):
            for j in range(4):M[idx[i]+axis][idx[j]+axis]=Maxis[i][j]
    return M


def _append_ba(M18,ba_margin:float):
    M=_zero(21)
    for i in range(18):
        for j in range(18):M[i][j]=M18[i][j]
    for i in range(18,21):M[i][i]=I(float(ba_margin))
    return M


def raw_blocks(domain_path:Path=HPF.DEFAULT_DOMAIN):
    path=Path(domain_path).resolve();h=H18FA.build();a=A21FA.build();dyn=DYNAMIC.build();proc=PROCESS.build()
    bad={'H18':H18FA.validate(h),'A21':A21FA.validate(a),'dynamic':DYNAMIC.validate(dyn),'process':PROCESS.validate(proc)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('strict-block prerequisites failed: '+repr(bad))
    penalty=float(h['delta_squared_completion_penalty']);ba=float(a['first_active_ba_margin_lower']);rows=[]
    for x in HPF._x_cover(dyn):
        M18=_h18_matrix(x,process=proc,dynamic=dyn,penalty_physical=penalty);ok18,p18=symmetric_positive_definite_ldlt(M18)
        M21=_append_ba(M18,ba);ok21,p21=symmetric_positive_definite_ldlt(M21)
        if not(ok18 and ok21):raise RuntimeError('materialized finite-angle strict block failed LDLT on x='+repr(x.as_list()))
        rows.append({'x':x,'H18':M18,'A21':M21,'H18_pivots':p18,'A21_pivots':p21})
    if not rows:raise RuntimeError('empty universal strict-block cover')
    return rows,h,a


def build(domain_path:Path=HPF.DEFAULT_DOMAIN):
    rows,h,a=raw_blocks(domain_path)
    h_p=min(float(p.lo) for r in rows for p in r['H18_pivots']);a_p=min(float(p.lo) for r in rows for p in r['A21_pivots']);reported=float(h['worst_full_H18_LDLT_pivot_lower'])
    if not math.isclose(h_p,reported,rel_tol=1e-12,abs_tol=1e-30):raise RuntimeError(f'materialized H18 pivot {h_p} differs from reported {reported}')
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','cell_count':len(rows),'H18_dimension':18,'A21_dimension':21,
      'finite_angle_H18_actual_interval_matrices_materialized':True,'A21_exact_first_active_direct_sum_matrices_materialized':True,
      'same_x_cover_as_canonical_H18_completion':True,'same_completion_penalty_as_finite_angle_H18':True,'P3_delta':HPF.USEFUL_GATE,
      'finite_angle_information_retention_lower':h['finite_angle_vector_information_retention_lower'],'H18_worst_materialized_LDLT_pivot_lower':h_p,
      'H18_reported_LDLT_pivot_lower':reported,'A21_worst_materialized_LDLT_pivot_lower':a_p,'A21_first_active_ba_margin_lower':a['first_active_ba_margin_lower'],
      'raw_block_API':'raw_blocks(domain_path)->[{x,H18,A21,H18_pivots,A21_pivots},...]','scalar_margin_substituted_for_matrix':False,
      'source_enumeration_used':False,'trajectory_replay_used':False,'declared_domain_shrunk':False,'endpoint_augmented_LDLT_closed_here':False,
      'every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False}


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    for k in ('finite_angle_H18_actual_interval_matrices_materialized','A21_exact_first_active_direct_sum_matrices_materialized','same_x_cover_as_canonical_H18_completion','same_completion_penalty_as_finite_angle_H18'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('scalar_margin_substituted_for_matrix','source_enumeration_used','trajectory_replay_used','declared_domain_shrunk','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    if int(d.get('cell_count',0))<=0:f.append('empty strict-block cover')
    for k in ('H18_worst_materialized_LDLT_pivot_lower','A21_worst_materialized_LDLT_pivot_lower','A21_first_active_ba_margin_lower'):
        if not(math.isfinite(float(d.get(k,math.nan))) and float(d[k])>0):f.append(k+' invalid')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=HPF.DEFAULT_DOMAIN);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'cells':d['cell_count'],'H18_pivot':d['H18_worst_materialized_LDLT_pivot_lower'],'A21_pivot':d['A21_worst_materialized_LDLT_pivot_lower'],'A21_ba':d['A21_first_active_ba_margin_lower'],'failures':f},sort_keys=True));return int(bool(f))

if __name__=='__main__':raise SystemExit(main())
