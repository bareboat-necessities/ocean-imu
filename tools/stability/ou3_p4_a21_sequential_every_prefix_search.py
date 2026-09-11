#!/usr/bin/env python3
"""All-family A21 every-literal-prefix sequential Schur search.

The canonical A21 post-prediction word has one new in-word physical/bias source
block only at its terminal following prediction.  Earlier literal prefixes
therefore carry the same persistent true-bias state/sector and event-local
binary32 inputs, but no new BRMM acceleration or [w,m_tau] driver block.  The
terminal prediction prefix is exactly the endpoint problem and delegates to the
all-family endpoint search.

For every BIAS0/BIAS1/BIAS2 family this module searches a finite Gamma_l and
source/noise coefficients against the actual shipping covariance information
at that prefix.  Point success is diagnostic only and never promotes P4.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_sub,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_a21_post_prediction_word as WORD
import ou3_p4_a21_source_indexed_prefix_transport as A21
import ou3_p4_a21_sequential_source_bias_fp_schur as SCHUR
import ou3_p4_a21_sequential_word_search as ENDPOINT
import ou3_p4_bias_family_joint_iss_supply as BIAS
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA=1
QUALIFICATION='OU3_P4_A21_ALL_BIAS_SEQUENTIAL_EVERY_LITERAL_PREFIX_SEARCH_V1'
P3_DELTA=1e-18
NE=21;NP=SCHUR.NP

def I(x):return Interval.point(float(x))
def zero(n):return [[I(0) for _ in range(n)] for _ in range(n)]
def embed_information(J):
    Q=zero(NP)
    for i in range(NE):
        for j in range(NE):Q[i][j]=J[i][j]
    return Q

def gain_margin(required,J0,Gamma):
    if not(math.isfinite(Gamma) and Gamma>0):raise ValueError('positive finite Gamma required')
    D=zero(NP)
    for i in range(NE):
        for j in range(NE):D[i][j]=I(Gamma)*J0[i][j]
    M=matrix_symmetric_hull(matrix_sub(D,required));ok,piv=symmetric_positive_definite_ldlt(M)
    return ok,[float(x.lo) for x in piv]

def backward_ordinary_prefix(word,prefix_len,family,*,gamma_n,beta_multiplier):
    if not(1<=prefix_len<len(word['events'])):raise ValueError('ordinary prefix must precede terminal prediction')
    events=word['events'][:prefix_len]
    if any(e['kind']=='prediction' for e in events):raise ValueError('ordinary prefix unexpectedly contains prediction')
    fp=FPISS.build();ff=FPISS.validate(fp);bs=BIAS.build();bf=BIAS.validate(bs)
    if ff or bf:raise RuntimeError(f'prefix prerequisites failed fp={ff} bias={bf}')
    beta=float(bs['family_supply'][family]['true_bias_norm_upper_mps2'])
    Qnext=embed_information(matrix_inverse_gauss_jordan(word['covariance_nodes'][prefix_len]));local=[]
    for rev,event in enumerate(reversed(events),start=1):
        delta=ENDPOINT._local_delta(event['kind'],fp)
        T=SCHUR.ordinary_transition(event['C'],delta);sectors=SCHUR.ordinary_sectors(len(T[0]),beta)
        try:
            K=SCHUR.local_master(Qnext,T,gamma_phys=0.0,gamma_bias=0.0,gamma_n=gamma_n,sectors=sectors,multipliers=(beta_multiplier,),prediction=False)
            Qprev,_,piv=SCHUR.schur_required_state_form(K)
        except ValueError as exc:
            return {'closed':False,'failure_class':'C','failure_stage':'local_Kuu_positive_definiteness','kind':event['kind'],'reverse_ordinal':rev,'reason':str(exc),'local_records':local}
        local.append({'kind':event['kind'],'Kuu_pivot_lower_min':min(piv)});Qnext=Qprev
    return {'closed':True,'required':Qnext,'J0':matrix_inverse_gauss_jordan(word['entrance_P']),'local_records':local,'true_bias_norm':beta}

def search_ordinary_prefix(word,prefix_len,family):
    vals=(1e-6,1e-3,1.0,1e3,1e6,1e9,1e12,1e15,1e18)
    Gammas=(1.0,1.25,1.5,2.0,3.0,5.0,10.0,30.0,100.0,1e3,1e4,1e6)
    attempts=0;best=None;first_local=None
    for gn in vals:
      for lbeta in vals:
        rec=backward_ordinary_prefix(word,prefix_len,family,gamma_n=gn,beta_multiplier=lbeta);attempts+=1
        if not rec['closed']:
            if first_local is None:first_local={k:v for k,v in rec.items() if k not in ('local_records','required')}
            continue
        for G in Gammas:
            ok,piv=gain_margin(rec['required'],rec['J0'],G);margin=min(piv) if piv else -math.inf
            row={'Gamma':G,'gamma_n':gn,'beta_multiplier':lbeta,'pivot_lower_min':margin,'LDLT_closed':ok,'local_Kuu_pivot_lower_min':min(x['Kuu_pivot_lower_min'] for x in rec['local_records']),'true_bias_norm':rec['true_bias_norm']}
            if best is None or margin>best['pivot_lower_min']:best=row
            if ok:return {'closed':True,'attempts':attempts,'candidate':row,'first_local_failure':first_local}
    return {'closed':False,'attempts':attempts,'candidate':best,'first_local_failure':first_local,'classification':'C'}

def audit_family(family):
    samples=A21.smoke_objects(family);word=WORD.build_word(samples,family);rows=[]
    for ell,event in enumerate(word['events'],start=1):
        if ell==len(word['events']):
            ep=ENDPOINT.search_family(samples,family)
            candidate=ep.get('candidate')
            r={'closed':bool(ep.get('candidate_found')),'attempts':ep.get('attempts',0),'candidate':candidate,'first_local_failure':ep.get('first_local_failure'),'endpoint_problem_reused':True}
        else:
            r=search_ordinary_prefix(word,ell,family);r['endpoint_problem_reused']=False
        rows.append({'prefix':ell,'kind':event['kind'],**r})
    closed=all(r['closed'] for r in rows);worst=next((r for r in rows if not r['closed']),None)
    if closed:
        def gain(r):
            c=r.get('candidate') or {}
            return float(c.get('Gamma',c.get('rho',0.0)))
        worst=max(rows,key=gain)
    return {'all_prefixes_have_finite_point_bound':closed,'prefix_count':len(rows),'worst_prefix':worst,'prefixes':rows}

def build():
    w=WORD.build();wf=WORD.validate(w);op=SCHUR.build();of=SCHUR.validate(op);ep=ENDPOINT.build();ef=ENDPOINT.validate(ep)
    if wf or of or ef:raise RuntimeError(f'A21 prefix-search prerequisites failed word={wf} op={of} endpoint={ef}')
    rows={f:audit_family(f) for f in BIAS.REQUIRED_BIAS_FAMILIES}
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'all_BIAS0_BIAS1_BIAS2_literal_prefixes_searched':True,'actual_prefix_covariance_information_targets_consumed':True,'preterminal_prefixes_do_not_double_charge_future_prediction_source':True,'persistent_true_bias_sector_retained_preterminal':True,'conditional_binary32_inputs_eliminated_event_by_event':True,'terminal_prefix_reuses_same_all_family_endpoint_problem':True,
      'family_prefix_search':rows,'point_prefix_bounds_are_not_source_uniform_certificates':True,'production_every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'freeze only successful prefix multiplier/gain topology; certify every prefix over dependency-preserving source cells, then combine finite Gamma/kappa values with coercivity and physical chart radii in the first-exit inequalities'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('all_BIAS0_BIAS1_BIAS2_literal_prefixes_searched','actual_prefix_covariance_information_targets_consumed','preterminal_prefixes_do_not_double_charge_future_prediction_source','persistent_true_bias_sector_retained_preterminal','conditional_binary32_inputs_eliminated_event_by_event','terminal_prefix_reuses_same_all_family_endpoint_problem','point_prefix_bounds_are_not_source_uniform_certificates'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    rows=d.get('family_prefix_search',{})
    if set(rows)!=set(BIAS.REQUIRED_BIAS_FAMILIES):f.append('family prefix table incomplete')
    for name,row in rows.items():
        if not isinstance(row.get('prefix_count'),int) or row['prefix_count']<=0:f.append(name+' prefix search did not execute')
        if len(row.get('prefixes',[]))!=row.get('prefix_count'):f.append(name+' prefix count mismatch')
    return f

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'families':{k:{'all_prefixes':v['all_prefixes_have_finite_point_bound'],'count':v['prefix_count'],'worst':v['worst_prefix']} for k,v in d['family_prefix_search'].items()},'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
