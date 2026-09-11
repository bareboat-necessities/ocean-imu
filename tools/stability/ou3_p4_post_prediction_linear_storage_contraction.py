#!/usr/bin/env python3
"""Exact consecutive-word linear storage contraction inherited from P3.

Canonical conditional P3 proves the full-matrix moving-information inequality at

    delta = 10^-18

for H18 and for the H18->A21/A21 continuation.  The P4 storage is now chosen at
the same post-prediction word boundaries as

    W_k = e_k^T P_k^{-1} e_k.

The boundary covariance is the identical shipping object inherited by the
successor word, so the metric transition penalty is mu=1.  Therefore the
homogeneous linear core satisfies exactly

    W_{k+1} <= rho_0 W_k,
    rho_0 = 1-delta = 999999999999999999 / 1000000000000000000 < 1.

The factor is retained as an exact integer ratio.  Computing `1.0-1e-18` in
binary64 produces 1.0 and is forbidden for theorem promotion.

This is the strict core around which nonlinear graph sectors and additive
physical/bias/binary32 supplies must be closed.  It is not itself P4.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

import ou3_brmm_riccati_metric_p3 as P3
import ou3_p4_information_storage_coercivity as STORAGE

SCHEMA=1
QUALIFICATION='OU3_P4_EXACT_POST_PREDICTION_LINEAR_STORAGE_CONTRACTION_V1'
DELTA=Fraction(1,10**18)
RHO=Fraction(10**18-1,10**18)

def build():
    p=P3.build();pf=P3.validate(p);s=STORAGE.build();sf=STORAGE.validate(s)
    if pf or sf:raise RuntimeError(f'linear storage prerequisites failed P3={pf} storage={sf}')
    h=bool(p['modes']['H18']['Omega_minus_delta_P_full_matrix_closed'])
    a=bool(p['modes']['A21']['Omega_minus_delta_P_full_matrix_closed'])
    exact=bool(p['P3_CONDITIONAL_BRMM_PASS'] and h and a and s['uniform_word_boundary_storage_coercivity_closed'] and s['consecutive_word_metric_comparison_mu']==1.0)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'P3_delta_exact_numerator':DELTA.numerator,'P3_delta_exact_denominator':DELTA.denominator,
      'rho0_exact_numerator':RHO.numerator,'rho0_exact_denominator':RHO.denominator,
      'rho0_decimal_string':'0.'+'9'*18,
      'rho0_strictly_between_zero_and_one':0<RHO<1,
      'binary64_one_minus_delta_may_be_used_for_strictness':False,
      'H18_full_matrix_P3_contraction_consumed':h,'A21_full_matrix_P3_contraction_consumed':a,
      'post_prediction_uniform_storage_coercivity_consumed':True,
      'consecutive_word_metric_transition_penalty_mu_exact':'1',
      'linear_consecutive_word_storage_contraction_closed':exact,
      'nonlinear_graph_sectors_absorbed_here':False,'physical_source_supply_absorbed_here':False,'bias_supply_absorbed_here':False,'binary32_supply_absorbed_here':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'start the augmented endpoint master from this exact strict rational core; add same-history nonlinear graph sectors and correlated BRMM/BIAS/binary32 supplies without converting the 1e-18 strict margin to binary64'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta_exact_numerator')!=1 or d.get('P3_delta_exact_denominator')!=10**18:f.append('P3 delta changed')
    if d.get('rho0_exact_numerator')!=10**18-1 or d.get('rho0_exact_denominator')!=10**18:f.append('rho0 exact ratio changed')
    for k in ('rho0_strictly_between_zero_and_one','H18_full_matrix_P3_contraction_consumed','A21_full_matrix_P3_contraction_consumed','post_prediction_uniform_storage_coercivity_consumed','linear_consecutive_word_storage_contraction_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('binary64_one_minus_delta_may_be_used_for_strictness','nonlinear_graph_sectors_absorbed_here','physical_source_supply_absorbed_here','bias_supply_absorbed_here','binary32_supply_absorbed_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('consecutive_word_metric_transition_penalty_mu_exact')!='1':f.append('metric penalty changed')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'rho':[d['rho0_exact_numerator'],d['rho0_exact_denominator']],'strict':d['linear_consecutive_word_storage_contraction_closed'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
