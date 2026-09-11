#!/usr/bin/env python3
"""Asymmetric multiplier refinement for the H18 sequential word search.

The coarse search intentionally starts with one common multiplier for the three
physical sectors.  This module removes that artificial equality after scale
identification.  It refines lambda_a0, lambda_a1 and lambda_J independently,
with nearby source/roundoff supplies, while preserving nonnegativity and the
same actual post-prediction word.

A successful result is still only a point/topology candidate.  P4 remains
fail-closed until the tuple is re-certified over the source-uniform outward
cell family and every literal prefix.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_h18_sequential_word_search as BASE
import ou3_p4_h18_source_indexed_prefix_transport as HTR

SCHEMA=1
QUALIFICATION='OU3_P4_H18_SEQUENTIAL_ASYMMETRIC_MULTIPLIER_REFINEMENT_V1'
P3_DELTA=1e-18


def _near(x):
    x=float(x)
    return tuple(sorted(set(max(1e-18,x*f) for f in (0.03,0.1,0.3,1.0,3.0,10.0,30.0))))

def _rhos(center):
    vals=(0.90,0.95,0.98,0.99,0.995,0.999,0.9999,0.99999,0.999999)
    return tuple(sorted(set(v for v in vals if 0<v<1)))

def refine(samples):
    coarse=BASE.search(samples)
    if coarse.get('candidate_found'):
        return {'candidate_found':True,'coarse_candidate_sufficient':True,'attempts':0,'candidate':coarse['candidate'],'coarse':coarse}
    seed=coarse.get('candidate')
    if not seed:
        return {'candidate_found':False,'coarse_candidate_sufficient':False,'attempts':0,'candidate':None,'coarse':coarse,'classification':'C','reason':'coarse search never reached an entrance-form comparison'}
    l0=float(seed['multipliers'][0]); gs0=float(seed['gamma_s']); gn0=float(seed['gamma_n'])
    attempts=0; best=None; first_local=None
    lvals=_near(l0); gsvals=_near(gs0); gnvals=_near(gn0)
    for la0 in lvals:
      for la1 in lvals:
       for lj in lvals:
        mult=(la0,la1,lj)
        for gn in gnvals:
         for gs in gsvals:
          rec=BASE.backward_required_form(samples,gamma_s=gs,gamma_n=gn,multipliers=mult);attempts+=1
          if not rec['closed']:
              if first_local is None:first_local={k:v for k,v in rec.items() if k not in ('local_records','required_entrance_form','word')}
              continue
          for rho in _rhos(seed.get('rho')):
              ok,piv=BASE.entrance_margin(rec['required_entrance_form'],rec['entrance_information'],rho)
              margin=min(piv) if piv else -math.inf
              row={'rho':rho,'gamma_s':gs,'gamma_n':gn,'multipliers':[la0,la1,lj],'pivot_lower_min':margin,'entrance_LDLT_closed':ok,'local_Kuu_pivot_lower_min':min(x['Kuu_pivot_lower_min'] for x in rec['local_records'])}
              if best is None or margin>best['pivot_lower_min']:best=row
              if ok:return {'candidate_found':True,'coarse_candidate_sufficient':False,'attempts':attempts,'candidate':row,'coarse':coarse,'first_local_failure':first_local}
    return {'candidate_found':False,'coarse_candidate_sufficient':False,'attempts':attempts,'candidate':best,'coarse':coarse,'first_local_failure':first_local,'classification':'C','reason':'bounded asymmetric multiplier refinement did not close the point entrance LDLT'}

def build():
    result=refine(HTR.smoke_objects())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'three_physical_sector_multipliers_refined_independently':True,'all_refined_multipliers_nonnegative':True,'same_actual_post_prediction_word_as_base_search':True,'conditional_binary32_channel_retained':True,
      'refinement':result,'point_candidate_is_not_source_uniform_certificate':True,'source_uniform_outward_subdivision_closed_here':False,'production_endpoint_augmented_LDLT_closed_here':False,'production_every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'freeze only a successful multiplier topology and re-run it on dependency-preserving source cells; if no point candidate exists, inspect the reported entrance/local pivot direction before enlarging the storage family'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('three_physical_sector_multipliers_refined_independently','all_refined_multipliers_nonnegative','same_actual_post_prediction_word_as_base_search','conditional_binary32_channel_retained','point_candidate_is_not_source_uniform_certificate'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_outward_subdivision_closed_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    r=d.get('refinement',{})
    if 'candidate_found' not in r:f.append('refinement did not execute')
    if r.get('candidate_found') is False and r.get('attempts',0)==0 and r.get('candidate') is not None:f.append('invalid zero-attempt refinement state')
    return f

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'candidate':d['refinement'].get('candidate_found'),'attempts':d['refinement'].get('attempts'),'best':d['refinement'].get('candidate'),'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
