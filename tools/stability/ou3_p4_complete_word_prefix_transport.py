#!/usr/bin/env python3
"""Exact every-prefix transport family for complete-BRMM P4.

The complete-word endpoint transport identity is prefix stable.  For every
literal event prefix ell, apply the *same* variation-of-constants construction
to events[0:ell], embeddings[0:ell+1], epsilon[0:ell+1].  This yields

    M_ell,
    d_ell = r_ell + E_ell eps_ell - M_ell E_0 eps_0
            + sum_{k<ell} M_{ell:k+1}(C_k-L_k)E_k eps_k.

No endpoint substitution is allowed: production P4 must build and certify this
object for every literal prefix.  Because the truncated sequence retains the
same event records, every due S=0 operation inside a prefix uses its actual
applied R_S and every later suffix map is exactly the suffix for that prefix.

This module provides the exact algebraic prefix-family constructor and an exact
rational regression.  It does not itself provide the source-uniform interval
event sequence or close any LDLT.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

import ou3_p4_complete_word_endpoint_transport as ENDPOINT

SCHEMA=1
QUALIFICATION='OU3_P4_COMPLETE_BRMM_EVERY_PREFIX_TRANSPORT_V1'

def prefix_decompositions(events,embeddings,eps_nodes):
    if not events:raise ValueError('nonempty event sequence required')
    if len(embeddings)!=len(events)+1 or len(eps_nodes)!=len(events)+1:
        raise ValueError('one embedding/epsilon node required per boundary')
    out=[]
    for ell in range(1,len(events)+1):
        d=ENDPOINT.endpoint_decomposition(events[:ell],embeddings[:ell+1],eps_nodes[:ell+1])
        out.append({'prefix_length':ell,'terminal_kind':events[ell-1].get('kind'),'transport':d})
    return out

def _zero_residual(v):return all(x==0 for x in v)

def _smoke():
    F=Fraction
    I2=[[F(1),F(0)],[F(0),F(1)]]
    E0=[[F(1)],[F(0)]]
    E1=[[F(1)],[F(0)]]
    E2=[[F(1)],[F(0)]]
    # Event 0 is prediction/source-like: C=L, so interior shift cancels.
    C0=[[F(9,10),F(1,20)],[F(0),F(19,20)]]
    # Event 1 is accelerometer-like: C-L nonzero and must survive with suffix.
    L1=I2
    C1=[[F(4,5),F(0)],[F(0),F(1)]]
    events=[
      {'C':C0,'L':C0,'rho':[F(1,100),F(-1,200)],'kind':'prediction'},
      {'C':C1,'L':L1,'rho':[F(1,300),F(1,400)],'kind':'accelerometer'},
    ]
    eps=[[F(1,10)],[F(1,12)],[F(1,14)]]
    rec=prefix_decompositions(events,[E0,E1,E2],eps)
    return rec

def build():
    e=ENDPOINT.build();ef=ENDPOINT.validate(e)
    if ef:raise RuntimeError('endpoint transport prerequisite failed: '+repr(ef))
    rec=_smoke()
    exact=all(_zero_residual(x['transport']['identity_residual']) for x in rec)
    first_no_interior=rec[0]['transport']['interior_event_indices']==[]
    second_acc_only=rec[1]['transport']['interior_event_indices']==[1]
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'whole_word_endpoint_transport_reused_without_new_algebra':True,
      'literal_prefix_truncation_semantics':True,'every_prefix_transport_family_constructor_available':True,
      'prefix_M_and_defect_decomposition_emitted':True,'endpoint_cannot_substitute_for_prefixes':True,
      'actual_applied_RS_suffix_semantics_inherited':bool(e['actual_RS_regularization_enters_every_applicable_suffix']),
      'prediction_source_hybrid_floor_cancellation_inherited':bool(e['prediction_source_hybrid_floor_interior_epsilon_terms_cancel_exactly']),
      'accelerometer_only_interior_shift_semantics_inherited':bool(e['accepted_accelerometer_is_only_interior_epsilon_event_class']),
      'exact_rational_every_prefix_identity_closed':exact,'exact_rational_prediction_prefix_has_no_interior_shift':first_no_interior,
      'exact_rational_accelerometer_prefix_retains_interior_shift':second_acc_only,
      'source_uniform_interval_event_sequence_materialized_here':False,'production_every_prefix_transport_enclosed_here':False,
      'production_every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'feed each estimator-owned interval event sequence into prefix_decompositions, then pull the structured chord/reset/source/fp graphs through each returned prefix transport and certify every prefix augmented LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('whole_word_endpoint_transport_reused_without_new_algebra','literal_prefix_truncation_semantics','every_prefix_transport_family_constructor_available','prefix_M_and_defect_decomposition_emitted','endpoint_cannot_substitute_for_prefixes','actual_applied_RS_suffix_semantics_inherited','prediction_source_hybrid_floor_cancellation_inherited','accelerometer_only_interior_shift_semantics_inherited','exact_rational_every_prefix_identity_closed','exact_rational_prediction_prefix_has_no_interior_shift','exact_rational_accelerometer_prefix_retains_interior_shift'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_interval_event_sequence_materialized_here','production_every_prefix_transport_enclosed_here','production_every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'prefix_transport':d['every_prefix_transport_family_constructor_available'],'exact':d['exact_rational_every_prefix_identity_closed'],'production':d['production_every_prefix_transport_enclosed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
