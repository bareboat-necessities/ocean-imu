#!/usr/bin/env python3
"""Structural interval-safe wrapper for complete-word/prefix Phi transport.

The canonical endpoint identity is algebraic.  Blind interval subtraction can
lose that algebra: for a non-point Interval X, ``X-X`` is not [0,0].  Therefore
an interval implementation must not rediscover known exact identities such as
C=L by subtracting two copies of the same enclosure.

Production event records may provide either:

* ``C_minus_L``: a directly constructed enclosure of the exact algebraic
  difference (preferred for Joseph/reset, where C-L=-G K H); or
* ``C_equals_L_exact=True``: a structural witness that the exact difference is
  zero (prediction, source rebase, hybrid lift, covariance floor/identity).

This module evaluates the same variation-of-constants formula as
``ou3_p4_complete_word_endpoint_transport`` while preserving those structural
identities.  It also creates exact zero/one scalar constants for Interval input
instead of forming x-x.  No theorem assumption is added: callers are responsible
for constructing the witnesses from the actual shipping operation.
"""
from __future__ import annotations

import argparse,json
from fractions import Fraction
from pathlib import Path
from typing import Any,Sequence

from ou3_interval import Interval
import ou3_p4_complete_word_endpoint_transport as END

SCHEMA=1
QUALIFICATION='OU3_P4_STRUCTURAL_INTERVAL_PREFIX_TRANSPORT_V1'
P3_DELTA=1.0e-18


def _shape(A):return END._shape(A)
def _exact_zero_like(x):return Interval.point(0.0) if isinstance(x,Interval) else x-x
def _exact_one_like(x):return Interval.point(1.0) if isinstance(x,Interval) else x-x+1
def _zeros(r,c,z):return [[z for _ in range(c)] for _ in range(r)]
def _zvec(n,z):return [z for _ in range(n)]
def _is_exact_zero_matrix(A,z):return all(x==z for row in A for x in row)
def _contains_zero(x):
    if isinstance(x,Interval):return x.lo<=0.0<=x.hi
    return x==0

def _difference_map(event,z):
    cr,cc=_shape(event['C'])
    if event.get('C_minus_L') is not None:
        D=event['C_minus_L']
        if _shape(D)!=(cr,cc):raise ValueError('C_minus_L shape mismatch')
        return D
    if event.get('C_equals_L_exact') is True:return _zeros(cr,cc,z)
    return END._msub(event['C'],event['L'])

def endpoint_decomposition(events,embeddings,eps_nodes):
    n=len(events)
    if not n:raise ValueError('nonempty event sequence required')
    if len(embeddings)!=n+1 or len(eps_nodes)!=n+1:raise ValueError('one E/epsilon per boundary required')
    if not eps_nodes[0]:raise ValueError('nonempty epsilon coordinate required')
    z=_exact_zero_like(eps_nodes[0][0]);one=_exact_one_like(eps_nodes[0][0])
    suffix=END.suffix_products(events,one,z)
    final_dim=_shape(events[-1]['C'])[0]
    direct=_zvec(final_dim,z);rword=_zvec(final_dim,z);interior=_zvec(final_dim,z)
    accel=[];indices=[]
    for k,event in enumerate(events):
        cr,cc=_shape(event['C']);lr,lc=_shape(event['L'])
        if (cr,cc)!=(lr,lc):raise ValueError('C/L shape mismatch')
        if len(event['rho'])!=cr:raise ValueError('rho dimension mismatch')
        incoming=END._mv(END._mm(event['L'],embeddings[k]),eps_nodes[k])
        outgoing=END._mv(embeddings[k+1],eps_nodes[k+1])
        xi=END._add(event['rho'],END._sub(outgoing,incoming))
        direct=END._add(direct,END._mv(suffix[k+1],xi))
        rword=END._add(rword,END._mv(suffix[k+1],event['rho']))
        D=_difference_map(event,z)
        weighted=END._mm(suffix[k+1],END._mm(D,embeddings[k]))
        if not _is_exact_zero_matrix(weighted,z):
            indices.append(k);interior=END._add(interior,END._mv(weighted,eps_nodes[k]))
            if event.get('kind')!='accelerometer':raise ValueError(f"non-accelerometer event {k} has nonzero structural interior epsilon transport")
            accel.append([[-x for x in row] for row in weighted])
    M=suffix[0]
    endpoint=END._sub(END._mv(embeddings[-1],eps_nodes[-1]),END._mv(END._mm(M,embeddings[0]),eps_nodes[0]))
    decomposed=END._add(rword,END._add(endpoint,interior));residual=END._sub(direct,decomposed)
    return {'M_word':M,'direct_defect':direct,'r_word':rword,'endpoint_epsilon_term':endpoint,
            'interior_epsilon_term':interior,'interior_event_indices':indices,
            'accelerometer_suffix_blocks':accel,'decomposed_defect':decomposed,'identity_residual':residual,
            'identity_residual_contains_zero':all(_contains_zero(x) for x in residual)}
def prefix_decompositions(events,embeddings,eps_nodes):
    if not events:raise ValueError('nonempty events required')
    return [{'prefix_length':ell,'terminal_kind':events[ell-1].get('kind'),'transport':endpoint_decomposition(events[:ell],embeddings[:ell+1],eps_nodes[:ell+1])} for ell in range(1,len(events)+1)]
def _smoke():
    # Non-point intervals are deliberate: naive C-L would be nonzero.
    x=Interval(.8,.9);z=Interval.point(0);o=Interval.point(1)
    C=[[x,z],[z,x]];E=[[[o],[z]],[[o],[z]]];eps=[[Interval(.1,.2)],[Interval(.15,.25)]]
    r=endpoint_decomposition([{'kind':'prediction','C':C,'L':C,'rho':[z,z],'C_equals_L_exact':True}],E,eps)
    rational=prefix_decompositions([{'kind':'prediction','C':[[Fraction(9,10)]],'L':[[Fraction(9,10)]],'rho':[Fraction(0)],'C_equals_L_exact':True}],[[[Fraction(1)]],[[Fraction(1)]]],[[Fraction(1,10)],[Fraction(1,12)]])[0]['transport']
    return {'interval_prediction_has_no_interior':r['interior_event_indices']==[],
            'interval_identity_residual_contains_zero':r['identity_residual_contains_zero'],
            'rational_identity_exact':rational['identity_residual']==[Fraction(0)]}
def build():
    e=END.build();ef=END.validate(e)
    if ef:raise RuntimeError('canonical endpoint prerequisite failed: '+repr(ef))
    s=_smoke();closed=all(s.values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'canonical_variation_of_constants_identity_retained':True,'interval_x_minus_x_not_used_for_structural_equalities':True,
      'explicit_C_minus_L_map_supported':True,'exact_C_equals_L_witness_supported':True,'Interval_zero_one_constants_exact':True,
      'structural_prefix_transport_closed':closed,'smoke':s,'new_physical_assumption_added':False,'packetwise_norm_budget_used':False,
      'production_source_uniform_prefix_LDLT_closed_here':False,'P4_PASS':False,'P5_MAY_START':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('canonical_variation_of_constants_identity_retained','interval_x_minus_x_not_used_for_structural_equalities','explicit_C_minus_L_map_supported','exact_C_equals_L_witness_supported','Interval_zero_one_constants_exact','structural_prefix_transport_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('new_physical_assumption_added','packetwise_norm_budget_used','production_source_uniform_prefix_LDLT_closed_here','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'structural':d['structural_prefix_transport_closed'],'smoke':d['smoke'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
