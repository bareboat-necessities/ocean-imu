#!/usr/bin/env python3
"""Every-literal-prefix sequential Schur search for the H18 P4 word.

For each literal prefix of the canonical post-prediction H18 word, construct the
actual target information metric from the shipping covariance node at that
prefix and backward-eliminate event-local binary32 inputs.  If the prefix
contains the following prediction, also eliminate the same-history physical
acceleration witness using the coupled a0/a1/J012 sectors.

This searches finite prefix gains Gamma_l and source/noise supplies.  A point
success is diagnostic only; production still requires one dependency-preserving
source-uniform cover and hard-domain first-exit inequalities.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_sub,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_h18_source_indexed_prefix_transport as HTR
import ou3_p4_h18_post_prediction_word as WORD
import ou3_p4_h18_sequential_source_fp_schur as SCHUR
import ou3_p4_h18_sequential_word_search as ENDPOINT
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA=1
QUALIFICATION='OU3_P4_H18_SEQUENTIAL_EVERY_LITERAL_PREFIX_SEARCH_V1'
P3_DELTA=1e-18
NX=SCHUR.NX;NP=SCHUR.NP

def I(x):return Interval.point(float(x))
def zero(n):return [[I(0) for _ in range(n)] for _ in range(n)]
def embed(J):
    Q=zero(NP)
    for i in range(NX):
        for j in range(NX):Q[i][j]=J[i][j]
    return Q

def by_token(samples):
    d={}
    for s in samples:
        for c in s.cells:d[c.source_token]=c
    return d

def backward_prefix(samples,w,prefix_len,*,gamma_s,gamma_n,multipliers):
    if not (1<=prefix_len<=len(w['events'])):raise ValueError('literal prefix outside word')
    fp=FPISS.build();ff=FPISS.validate(fp)
    if ff:raise RuntimeError('binary32 prerequisite open')
    cells=by_token(samples);Qnext=embed(matrix_inverse_gauss_jordan(w['covariance_nodes'][prefix_len]));local=[]
    events=w['events'][:prefix_len]
    for rev,event in enumerate(reversed(events),start=1):
        kind=event['kind'];delta=ENDPOINT._local_delta(kind,fp)
        if kind=='prediction':
            cell=cells.get(event['token'])
            if cell is None:raise RuntimeError('prediction token detached')
            h=float(cell.dt_s.hi)
            T=SCHUR.prediction_transition(event['C'],cell.tau_applied_s,cell.dt_s,delta)
            named=SCHUR.prediction_source_sectors(len(T[0]),h=h);sectors=[P for _,P in named];pred=True
        else:
            T=SCHUR.ordinary_transition(event['C'],delta);sectors=[];pred=False
        try:
            K=SCHUR.local_master(Qnext,T,gamma_s=gamma_s,gamma_n=gamma_n,sectors=sectors,multipliers=multipliers if pred else (),prediction=pred)
            Qprev,_,piv=SCHUR.schur_required_state_form(K)
        except ValueError as exc:
            return {'closed':False,'failure_class':'C','failure_stage':'local_Kuu_positive_definiteness','kind':kind,'reverse_ordinal':rev,'reason':str(exc),'local':local}
        local.append({'kind':kind,'Kuu_pivot_lower_min':min(piv)});Qnext=Qprev
    return {'closed':True,'required':Qnext,'J0':matrix_inverse_gauss_jordan(w['entrance_P']),'local':local,'contains_prediction':any(e['kind']=='prediction' for e in events)}

def gain_margin(required,J0,Gamma):
    if not(math.isfinite(Gamma) and Gamma>0):raise ValueError('positive finite Gamma required')
    desired=zero(NP)
    for i in range(NX):
        for j in range(NX):desired[i][j]=I(Gamma)*J0[i][j]
    M=matrix_symmetric_hull(matrix_sub(desired,required));ok,p=symmetric_positive_definite_ldlt(M)
    return ok,[float(x.lo) for x in p]

def search_one(samples,w,prefix_len):
    Gammas=(1.0,1.25,1.5,2.0,3.0,5.0,10.0,30.0,100.0,1e3,1e4,1e6)
    vals=tuple(10.0**p for p in range(-6,25,3));attempts=0;best=None;first_local=None
    contains_prediction=any(e['kind']=='prediction' for e in w['events'][:prefix_len])
    lambda_sets=((0.0,0.0,0.0),) if not contains_prediction else tuple((v,v,v) for v in vals)
    gs_vals=(0.0,) if not contains_prediction else vals
    for mult in lambda_sets:
      for gn in vals:
       for gs in gs_vals:
        rec=backward_prefix(samples,w,prefix_len,gamma_s=gs,gamma_n=gn,multipliers=mult);attempts+=1
        if not rec['closed']:
            if first_local is None:first_local={k:v for k,v in rec.items() if k!='local'}
            continue
        for G in Gammas:
            ok,p=gain_margin(rec['required'],rec['J0'],G);margin=min(p) if p else -math.inf
            row={'Gamma':G,'gamma_s':gs,'gamma_n':gn,'multipliers':list(mult),'pivot_lower_min':margin,'LDLT_closed':ok,'local_Kuu_pivot_lower_min':min(x['Kuu_pivot_lower_min'] for x in rec['local'])}
            if best is None or margin>best['pivot_lower_min']:best=row
            if ok:return {'closed':True,'attempts':attempts,'candidate':row,'contains_prediction':contains_prediction,'first_local_failure':first_local}
    return {'closed':False,'attempts':attempts,'candidate':best,'contains_prediction':contains_prediction,'first_local_failure':first_local,'classification':'C'}

def audit(samples):
    w=WORD.build_word(samples);rows=[]
    for ell,event in enumerate(w['events'],start=1):
        r=search_one(samples,w,ell);rows.append({'prefix':ell,'kind':event['kind'],**r})
    closed=all(r['closed'] for r in rows);worst=None
    if closed:
        worst=max(rows,key=lambda r:r['candidate']['Gamma'])
    else:
        worst=next(r for r in rows if not r['closed'])
    return {'all_prefixes_have_finite_point_bound':closed,'prefix_count':len(rows),'worst_prefix':worst,'prefixes':rows}

def build():
    base=WORD.build();bf=WORD.validate(base);sch=SCHUR.build();sf=SCHUR.validate(sch);fp=FPISS.build();ff=FPISS.validate(fp)
    if bf or sf or ff:raise RuntimeError(f'prefix search prerequisites failed word={bf} Schur={sf} fp={ff}')
    r=audit(HTR.smoke_objects())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'actual_literal_prefix_order_consumed':True,'actual_prefix_covariance_information_targets_consumed':True,'binary32_inputs_eliminated_event_by_event':True,'physical_source_sector_added_only_when_prefix_contains_prediction':True,'finite_Gamma_search_executed_for_every_literal_prefix':True,
      'diagnostic':r,'point_prefix_bounds_are_not_source_uniform_certificate':True,'production_every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use successful point prefix candidates only as multiplier/gain topology; certify the same inequalities over all dependency-preserving source cells and then compare physical prefix excursions against actual nonlinear chart/working-domain targets'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('actual_literal_prefix_order_consumed','actual_prefix_covariance_information_targets_consumed','binary32_inputs_eliminated_event_by_event','physical_source_sector_added_only_when_prefix_contains_prediction','finite_Gamma_search_executed_for_every_literal_prefix','point_prefix_bounds_are_not_source_uniform_certificate'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    r=d.get('diagnostic',{})
    if not isinstance(r.get('prefix_count'),int) or r['prefix_count']<=0:f.append('prefix search did not execute')
    if len(r.get('prefixes',[]))!=r.get('prefix_count'):f.append('prefix result count mismatch')
    return f

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'all_prefixes':d['diagnostic']['all_prefixes_have_finite_point_bound'],'count':d['diagnostic']['prefix_count'],'worst':d['diagnostic']['worst_prefix'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
