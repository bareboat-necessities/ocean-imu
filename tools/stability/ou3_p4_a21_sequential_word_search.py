#!/usr/bin/env python3
"""All-family A21 post-prediction sequential Schur search.

For each BIAS0/BIAS1/BIAS2 canonical A21 word, backward-eliminate the actual
literal-event binary32 inputs and, at the following prediction, the coupled
physical acceleration witness plus the family-specific shared [w,m_tau] bias
supply.  The persistent proof coordinate is the 29-state coordinate from
``ou3_p4_a21_sequential_source_bias_fp_schur``.

The endpoint storage is the shipping 21-state information metric embedded in
the joint coordinate.  The proof-only b_true/a_phys/source-scale coordinates
carry no artificial endpoint storage.  The true-bias hard family ball is
retained as a quadratic sector at every local event.

This is a point/topology search only.  Success does not promote P4; the tuple
must still be re-certified over the dependency-preserving source-uniform cells
and every literal prefix.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_sub,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_a21_post_prediction_word as WORD
import ou3_p4_a21_source_indexed_prefix_transport as A21
import ou3_p4_a21_sequential_source_bias_fp_schur as SCHUR
import ou3_p4_bias_family_joint_iss_supply as BIAS
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA=1
QUALIFICATION='OU3_P4_A21_ALL_BIAS_SEQUENTIAL_WORD_SEARCH_V1'
P3_DELTA=1e-18
NE=21;NP=SCHUR.NP

def I(x):return Interval.point(float(x))
def zero(n):return [[I(0) for _ in range(n)] for _ in range(n)]
def embed_information(J):
    if len(J)!=NE or any(len(r)!=NE for r in J):raise ValueError('A21 information matrix required')
    Q=zero(NP)
    for i in range(NE):
        for j in range(NE):Q[i][j]=J[i][j]
    return Q

def _local_delta(kind,fp):
    m=fp['modes']['A21'];dp=float(m['time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction']);dm=float(m['explicit_measurement_reset_state_roundoff_norm_upper_per_event'])
    if kind=='prediction':return dp
    if kind=='source_coordinate_rebase':return 0.0
    if kind=='aw_floor':return max(dp,dm)
    return dm

def _prediction_cell(samples,token):
    # A21 transport prediction token is the synchronized SourceCoverCell token.
    for s in samples:
        for c in s.cells:
            if c.source_token==token:return c
    raise RuntimeError('prediction token detached from A21 SourceCoverCell lineage')

def backward_required_form(samples,family,*,gamma_phys,gamma_bias,gamma_n,multipliers):
    if family not in BIAS.REQUIRED_BIAS_FAMILIES:raise ValueError('family must be BIAS0/BIAS1/BIAS2')
    if gamma_phys<0 or gamma_bias<0 or gamma_n<=0 or any(x<0 for x in multipliers):raise ValueError('invalid nonnegative supply/multiplier search value')
    word=WORD.build_word(samples,family);fp=FPISS.build();ff=FPISS.validate(fp);bs=BIAS.build();bf=BIAS.validate(bs)
    if ff or bf:raise RuntimeError(f'A21 search prerequisites failed fp={ff} bias={bf}')
    row=bs['family_supply'][family];supply_norm=float(row['joint_supply_norm_upper_per_prediction_mps2']);beta_norm=float(row['true_bias_norm_upper_mps2'])
    Qnext=embed_information(matrix_inverse_gauss_jordan(word['terminal_P']));local=[]
    for reverse_ordinal,event in enumerate(reversed(word['events']),start=1):
        kind=event['kind'];delta=_local_delta(kind,fp)
        if kind=='prediction':
            cell=_prediction_cell(samples,event['token']);h=float(cell.dt_s.hi)
            T=SCHUR.prediction_transition(event['C'],event['prediction_supply_block'],cell.tau_applied_s,cell.dt_s,delta)
            sectors=SCHUR.prediction_sectors(len(T[0]),h,supply_norm,beta_norm);pred=True
            if len(multipliers)!=len(sectors):raise ValueError('prediction multiplier count mismatch')
            lm=multipliers
        else:
            T=SCHUR.ordinary_transition(event['C'],delta);sectors=SCHUR.ordinary_sectors(len(T[0]),beta_norm);pred=False
            # one beta-sector multiplier for ordinary events, use the last entry
            lm=(multipliers[-1],)
        try:
            K=SCHUR.local_master(Qnext,T,gamma_phys=gamma_phys,gamma_bias=gamma_bias,gamma_n=gamma_n,sectors=sectors,multipliers=lm,prediction=pred)
            Qprev,_,piv=SCHUR.schur_required_state_form(K)
        except ValueError as exc:
            return {'closed':False,'failure_class':'C','failure_stage':'local_Kuu_positive_definiteness','kind':kind,'reverse_ordinal':reverse_ordinal,'reason':str(exc),'local_records':local}
        local.append({'kind':kind,'reverse_ordinal':reverse_ordinal,'Kuu_pivot_lower_min':min(piv)});Qnext=Qprev
    return {'closed':True,'required_entrance_form':Qnext,'entrance_information':matrix_inverse_gauss_jordan(word['entrance_P']),'word':word,'local_records':local,'family_supply_norm':supply_norm,'true_bias_norm':beta_norm}

def entrance_margin(required,J0,rho):
    if not(0<rho<1):raise ValueError('rho must be strictly below one')
    D=zero(NP)
    for i in range(NE):
        for j in range(NE):D[i][j]=I(rho)*J0[i][j]
    M=matrix_symmetric_hull(matrix_sub(D,required));ok,piv=symmetric_positive_definite_ldlt(M)
    return ok,[float(x.lo) for x in piv]

def grid():
    vals=tuple(10.0**p for p in range(-6,25,3));rhos=(.9,.95,.98,.99,.995,.999,.9999)
    return vals,rhos

def search_family(samples,family):
    vals,rhos=grid();attempts=0;best=None;first_local=None
    # grouped topology: three physical witness multipliers share one scale;
    # bias supply and beta-cap multipliers are independent.
    for lp in vals:
      for lbs in vals:
       for lbeta in vals:
        mult=(lp,lp,lp,lbs,lbeta)
        for gn in vals:
         for gp in vals:
          for gb in vals:
           rec=backward_required_form(samples,family,gamma_phys=gp,gamma_bias=gb,gamma_n=gn,multipliers=mult);attempts+=1
           if not rec['closed']:
             if first_local is None:first_local={k:v for k,v in rec.items() if k not in ('local_records','required_entrance_form','word')}
             continue
           for rho in rhos:
             ok,piv=entrance_margin(rec['required_entrance_form'],rec['entrance_information'],rho);margin=min(piv) if piv else -math.inf
             row={'rho':rho,'gamma_phys':gp,'gamma_bias':gb,'gamma_n':gn,'multipliers':list(mult),'pivot_lower_min':margin,'entrance_LDLT_closed':ok,'local_Kuu_pivot_lower_min':min(x['Kuu_pivot_lower_min'] for x in rec['local_records']),'family_supply_norm':rec['family_supply_norm'],'true_bias_norm':rec['true_bias_norm']}
             if best is None or margin>best['pivot_lower_min']:best=row
             if ok:return {'candidate_found':True,'attempts':attempts,'candidate':row,'first_local_failure':first_local,'event_count':len(rec['word']['events'])}
    return {'candidate_found':False,'attempts':attempts,'candidate':best,'first_local_failure':first_local,'classification':'C'}

def build():
    op=SCHUR.build();of=SCHUR.validate(op);word=WORD.build();wf=WORD.validate(word);bias=BIAS.build();bf=BIAS.validate(bias);fp=FPISS.build();ff=FPISS.validate(fp)
    if of or wf or bf or ff:raise RuntimeError(f'A21 word-search prerequisites failed op={of} word={wf} bias={bf} fp={ff}')
    rows={f:search_family(A21.smoke_objects(f),f) for f in BIAS.REQUIRED_BIAS_FAMILIES}
    successes=[f for f,r in rows.items() if r['candidate_found']]
    worst=None
    for f,r in rows.items():
        c=r.get('candidate')
        if c is None:continue
        if worst is None or c['pivot_lower_min']<worst['candidate']['pivot_lower_min']:worst={'family':f,'candidate':c}
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'all_three_bias_families_searched_separately':True,'same_shared_w_mtau_supply_used_in_search':True,'persistent_true_bias_sector_applied_at_every_event':True,'physical_and_bias_source_scales_kept_distinct':True,'conditional_binary32_event_channel_consumed':True,'strict_rho_below_one_search_executed':True,
      'family_search':rows,'families_with_point_candidate':successes,'worst_point_candidate':worst,'point_candidates_are_not_source_uniform_certificates':True,'source_uniform_outward_subdivision_closed_here':False,'production_endpoint_augmented_LDLT_closed_here':False,'production_every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use each family point result only to choose multiplier/storage topology; then certify the same endpoint and every-prefix inequalities over all dependency-preserving source cells, preserving projection branch and family ancestry'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('all_three_bias_families_searched_separately','same_shared_w_mtau_supply_used_in_search','persistent_true_bias_sector_applied_at_every_event','physical_and_bias_source_scales_kept_distinct','conditional_binary32_event_channel_consumed','strict_rho_below_one_search_executed','point_candidates_are_not_source_uniform_certificates'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_outward_subdivision_closed_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    rows=d.get('family_search',{})
    if set(rows)!=set(BIAS.REQUIRED_BIAS_FAMILIES):f.append('family search table incomplete')
    for name,row in rows.items():
        if not isinstance(row.get('attempts'),int) or row['attempts']<=0:f.append(name+' search did not execute')
    return f

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'successes':d['families_with_point_candidate'],'worst':d['worst_point_candidate'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
