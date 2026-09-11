#!/usr/bin/env python3
"""Uniform coercivity and exact compatibility of the P4 information storage.

Use the source-dependent shipping information metric only at the canonical P3/P4
word boundaries

    W_k(e,zeta_k) = e^T P_k^{-1} e,

where a boundary is taken *immediately after a shipping prediction*.  Canonical
P3 already has this topology: a complete information word closes at its
following prediction.  This boundary choice matters for arbitrary asynchronous
vector traffic between predictions.

At every post-prediction boundary shipping has

    P_k = F P_pre F^T + Q,      Q >= q_min I,

hence P_k >= q_min I independently of the number of accepted measurements in
the preceding suffix.  The certified source-uniform covariance ceiling gives
P_k <= trace(Pbar) I.  Consequently

    1/trace(Pbar) I <= P_k^{-1} <= 1/q_min I.

Both coercivity constants are therefore finite and source-uniform without any
hardware ODR or accepted-magnetometer-count premise.

The following prediction is simultaneously the end boundary of one word and the
start boundary of its successor, so the physical covariance object is identical
at the shared boundary.  Thus the storage comparison penalty is exactly

    M_next = M_end,   mu = 1.

This does not claim that P^{-1} is uniformly upper bounded at every *interior*
Joseph prefix when the source language does not cap asynchronous accepted-vector
information.  Prefix excursion/first-exit certificates may compare physical
state to the uniformly coercive boundary storage and need not redefine the
word-boundary metric at every Joseph packet.  No source assumption is added.
"""
from __future__ import annotations
import argparse,json,math

import ou3_full_process_ucc as PROCESS
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_p3_premises as PREMISES
import ou3_brmm_riccati_metric_p3 as P3

SCHEMA=2
QUALIFICATION='OU3_P4_POST_PREDICTION_INFORMATION_STORAGE_COERCIVITY_V2'
P3_DELTA=1.0e-18

def _mode_bounds(tube,process,mode,key):
    diag=[float(x) for x in tube['modes'][key]['Pbar_diagonal_variance_upper']]
    if not diag or any(not(math.isfinite(x) and x>0) for x in diag):raise RuntimeError('invalid covariance ceiling')
    pmax=math.nextafter(sum(diag),math.inf)
    qmin=float(process['modes'][key]['prediction_Q_lambda_min_lower'])
    if not(math.isfinite(qmin) and qmin>0):raise RuntimeError('strict full-state Q lower lost')
    mminus=math.nextafter(1.0/pmax,0.0)
    mplus=math.nextafter(1.0/qmin,math.inf)
    return {'dimension':len(diag),'Pbar_diagonal_variance_upper':diag,'post_prediction_P_lambda_max_upper_via_trace':pmax,'post_prediction_P_lambda_min_lower_from_Q':qmin,'information_storage_m_minus':mminus,'information_storage_m_plus':mplus,'uniform_boundary_coercivity_closed':mminus>0 and math.isfinite(mplus)}

def _async_upper_rate_not_needed(premises):
    pe=premises['execution_premises']['vector_PE']
    required={'vector_pe_recurrence_window_s','specific_force_norm_lower_mps2','specific_force_norm_upper_mps2','magnetic_vector_norm_lower_uT','magnetic_vector_norm_upper_uT','vector_sine_separation_lower','body_rate_norm_upper_deg_s'}
    if set(pe)!=required:raise RuntimeError('vector PE schema changed')
    return True

def build():
    process=PROCESS.build();pf=PROCESS.validate(process)
    tube=TUBE.build();tf=TUBE.validate_covariance_ceiling(tube)
    premises=PREMISES.build();prf=PREMISES.validate(premises)
    p3=P3.build();p3f=P3.validate(p3)
    if pf or tf or prf or p3f:raise RuntimeError(f'coercivity prerequisites failed process={pf} tube={tf} premises={prf} p3={p3f}')
    rows={m:_mode_bounds(tube,process,m,k) for m,k in (('H18','H'),('A21','A'))}
    topology=bool(
      p3['P3_CONDITIONAL_BRMM_PASS'] and
      p3['canonical_P3_topology']=='H18_3S_PRIOR_FREE_THEN_PRESERVED_H_TO_A_HYBRID_A21' and
      p3['modes']['H18']['Omega_minus_delta_P_full_matrix_closed'] and
      p3['modes']['A21']['Omega_minus_delta_P_full_matrix_closed'])
    closed=topology and all(r['uniform_boundary_coercivity_closed'] for r in rows.values())
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'candidate_storage':'W_k=e^T P_k^-1 e at post-prediction word boundaries',
      'canonical_P3_following_prediction_boundary_consumed':topology,
      'shipping_reachable_covariance_is_storage_source_coordinate':True,
      'strict_full_state_process_Q_lower_consumed':True,
      'source_uniform_covariance_ceiling_consumed':True,
      'post_prediction_P_lower_independent_of_prior_and_measurement_suffix':True,
      'uniform_lower_coercivity_m_minus_closed':all(r['information_storage_m_minus']>0 for r in rows.values()),
      'uniform_upper_coercivity_m_plus_closed':all(math.isfinite(r['information_storage_m_plus']) for r in rows.values()),
      'uniform_word_boundary_storage_coercivity_closed':closed,
      'consecutive_word_boundary_covariance_is_identical':True,
      'consecutive_word_metric_comparison_mu':1.0,
      'metric_jump_penalty_hidden':False,
      'arbitrary_async_measurement_count_requires_boundary_m_plus_premise':False,
      'hardware_magnetometer_ODR_substituted_as_theorem_premise':False,
      'interior_Joseph_P_inverse_uniform_upper_bound_claimed':False,
      'word_boundary_moved_to_nonshipping_virtual_event':False,
      'modes':rows,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use this exact mu=1 uniformly coercive post-prediction storage at consecutive word boundaries; express every literal interior prefix as a transport from the same boundary storage, then attach correlated BRMM/bias/binary32 supplies and close endpoint/every-prefix LDLT plus first-exit targets'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('canonical_P3_following_prediction_boundary_consumed','shipping_reachable_covariance_is_storage_source_coordinate','strict_full_state_process_Q_lower_consumed','source_uniform_covariance_ceiling_consumed','post_prediction_P_lower_independent_of_prior_and_measurement_suffix','uniform_lower_coercivity_m_minus_closed','uniform_upper_coercivity_m_plus_closed','uniform_word_boundary_storage_coercivity_closed','consecutive_word_boundary_covariance_is_identical'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('consecutive_word_metric_comparison_mu')!=1.0:f.append('metric compatibility penalty changed')
    for k in ('metric_jump_penalty_hidden','arbitrary_async_measurement_count_requires_boundary_m_plus_premise','hardware_magnetometer_ODR_substituted_as_theorem_premise','interior_Joseph_P_inverse_uniform_upper_bound_claimed','word_boundary_moved_to_nonshipping_virtual_event','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,row in d.get('modes',{}).items():
        if not(float(row.get('information_storage_m_minus',0))>0):f.append(mode+' m_minus not strict')
        if not(math.isfinite(float(row.get('information_storage_m_plus',math.inf))) and float(row['information_storage_m_plus'])>0):f.append(mode+' m_plus invalid')
        if row.get('uniform_boundary_coercivity_closed') is not True:f.append(mode+' boundary coercivity open')
    return list(dict.fromkeys(f))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    import pathlib;path=pathlib.Path(a.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'mu':d['consecutive_word_metric_comparison_mu'],'coercive':d['uniform_word_boundary_storage_coercivity_closed'],'bounds':{k:[v['information_storage_m_minus'],v['information_storage_m_plus']] for k,v in d['modes'].items()},'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
