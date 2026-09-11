#!/usr/bin/env python3
"""Coercivity and consecutive-word compatibility of reachable P^{-1} storage.

A natural P4 source-dependent storage is the shipping information metric

    W(e,zeta) = e^T P(zeta)^{-1} e.

It has an important advantage over a separately fitted Lyapunov family: the
shipping covariance at the end of one word is exactly the covariance inherited
by the next word.  Therefore the boundary metric is identical and the
inter-word comparison penalty is exactly mu=1; there is no hidden metric jump.

Uniform lower coercivity of P^{-1} follows from the already-certified covariance
ceiling P <= Pbar in the weak but sufficient form

    P <= trace(Pbar) I  =>  P^{-1} >= trace(Pbar)^{-1} I.

Uniform upper coercivity needs P >= p_min I.  The full shipping process UCC gives
Q >= q_min I at every IMU prediction, hence P^- >= q_min I independently of the
prior.  For exact linear Joseph updates,

    (P^+)^{-1} = (P^-)^{-1} + H^T R^{-1} H,

so after any finite accepted measurement suffix

    P^+ >= [q_min^{-1} + sum ||H_i||_2^2/lambda_min(R_i)]^{-1} I.

Thus a finite source-uniform upper bound on the *cumulative accepted measurement
information between consecutive predictions* is sufficient for uniform upper
coercivity.  The current Normal-Live/P3 asynchronous magnetic premise provides
recurring accepted packets and permits outages/rejections, but does not declare
an upper accepted-information rate or packet count.  The typed theorem source
also stores a tuple of post-IMU magnetic events without a source-uniform length
cap.  Therefore this P^{-1} storage route cannot yet claim a finite global m_+
from the declared premises alone.

This is not an instability result.  It is a missing execution/source
qualification (classification E) for this storage construction, or alternatively
a storage-architecture obligation (classification G) if a uniformly coercive
metric not requiring this information-rate premise is constructed.  The weakest
additional premise is a finite cumulative information bound, not a hardware ODR.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_p3_premises as PREMISES
import ou3_brmm_complete_window_execution_kernel as KERNEL

SCHEMA=1
QUALIFICATION='OU3_P4_REACHABLE_INFORMATION_STORAGE_COERCIVITY_V1'
P3_DELTA=1.0e-18

def _mode_ceiling(tube,mode):
    key='H' if mode=='H18' else 'A'
    diag=[float(x) for x in tube['modes'][key]['Pbar_diagonal_variance_upper']]
    if not diag or any(not(math.isfinite(x) and x>0) for x in diag):raise RuntimeError('invalid covariance ceiling')
    tr=math.nextafter(sum(diag),math.inf)
    return diag,tr,math.nextafter(1.0/tr,0.0)

def _declared_async_upper_rate_absent(premises):
    e=premises['execution_premises'];pe=e['vector_PE']
    # These are lower/recurrent observability and geometry assumptions.  None is
    # an accepted-packet-count or cumulative-information upper bound.
    required={'vector_pe_recurrence_window_s','specific_force_norm_lower_mps2','specific_force_norm_upper_mps2','magnetic_vector_norm_lower_uT','magnetic_vector_norm_upper_uT','vector_sine_separation_lower','body_rate_norm_upper_deg_s'}
    if set(pe)!=required:raise RuntimeError('vector PE schema changed; re-audit async information-rate conclusion')
    return True

def posterior_lower_from_information_budget(q_min:float,information_budget:float)->float:
    if not(math.isfinite(q_min) and q_min>0):raise ValueError('strict q_min required')
    if not(math.isfinite(information_budget) and information_budget>=0):raise ValueError('finite nonnegative information budget required')
    denom=math.nextafter(1.0/q_min+information_budget,math.inf)
    return math.nextafter(1.0/denom,0.0)

def build():
    process=PROCESS.build();pf=PROCESS.validate(process)
    tube=TUBE.build();tf=TUBE.validate_covariance_ceiling(tube)
    premises=PREMISES.build();prf=PREMISES.validate(premises)
    if pf or tf or prf:raise RuntimeError(f'coercivity prerequisites failed process={pf} tube={tf} premises={prf}')
    rows={}
    for mode,key in (('H18','H'),('A21','A')):
        diag,tr,mminus=_mode_ceiling(tube,mode)
        q=float(process['modes'][key]['prediction_Q_lambda_min_lower'])
        example=posterior_lower_from_information_budget(q,0.0)
        rows[mode]={
          'dimension':len(diag),'Pbar_diagonal_variance_upper':diag,
          'P_lambda_max_upper_via_trace':tr,
          'information_storage_coercivity_lower_m_minus':mminus,
          'prediction_Q_lambda_min_lower':q,
          'posterior_P_lower_if_information_budget_zero':example,
          'finite_m_plus_formula':'m_plus <= 1/q_min + I_meas,max per prediction interval',
          'finite_m_plus_closed_under_current_declared_premises':False}
    no_upper=_declared_async_upper_rate_absent(premises)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'candidate_storage':'W=e^T P(zeta)^-1 e','shipping_reachable_covariance_is_storage_source_coordinate':True,
      'consecutive_word_boundary_covariance_is_identical':True,
      'consecutive_word_metric_comparison_mu':1.0,
      'metric_jump_penalty_hidden':False,
      'uniform_lower_coercivity_m_minus_closed':all(r['information_storage_coercivity_lower_m_minus']>0 for r in rows.values()),
      'strict_full_state_process_Q_lower_consumed':True,
      'prediction_erases_need_for_prior_P_lower_bound':True,
      'Joseph_information_identity_consumed':True,
      'finite_cumulative_measurement_information_bound_is_sufficient_for_m_plus':True,
      'declared_async_magnetic_PE_has_no_uniform_accepted_information_rate_upper':no_upper,
      'hardware_magnetometer_ODR_substituted_as_theorem_premise':False,
      'uniform_upper_coercivity_m_plus_closed':False,
      'uniform_storage_coercivity_closed':False,
      'weakest_missing_premise':'finite source-uniform upper bound on sum_i ||H_i||_2^2/lambda_min(R_i) between consecutive IMU predictions (or an equivalent cumulative accepted-vector-information bound)',
      'failure_classification_if_using_information_storage':'E unless an alternative uniformly coercive compatible storage removes this premise, then G is the remaining architecture obligation',
      'modes':rows,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'either qualify the weakest cumulative accepted-information upper bound from the declared execution/source contract and obtain m_plus, or construct a different source-dependent storage with certified uniform coercivity and successor compatibility; do not assume hardware ODR'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('shipping_reachable_covariance_is_storage_source_coordinate','consecutive_word_boundary_covariance_is_identical','uniform_lower_coercivity_m_minus_closed','strict_full_state_process_Q_lower_consumed','prediction_erases_need_for_prior_P_lower_bound','Joseph_information_identity_consumed','finite_cumulative_measurement_information_bound_is_sufficient_for_m_plus','declared_async_magnetic_PE_has_no_uniform_accepted_information_rate_upper'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('consecutive_word_metric_comparison_mu')!=1.0:f.append('metric compatibility penalty is not exactly one')
    for k in ('metric_jump_penalty_hidden','hardware_magnetometer_ODR_substituted_as_theorem_premise','uniform_upper_coercivity_m_plus_closed','uniform_storage_coercivity_closed','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,row in d.get('modes',{}).items():
        if not(float(row.get('information_storage_coercivity_lower_m_minus',0))>0):f.append(mode+' m_minus not strict')
        if not(float(row.get('prediction_Q_lambda_min_lower',0))>0):f.append(mode+' q_min not strict')
        if row.get('finite_m_plus_closed_under_current_declared_premises') is not False:f.append(mode+' m_plus unexpectedly closed')
    return list(dict.fromkeys(f))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'mu':d['consecutive_word_metric_comparison_mu'],'m_minus':{k:v['information_storage_coercivity_lower_m_minus'] for k,v in d['modes'].items()},'m_plus_closed':d['uniform_upper_coercivity_m_plus_closed'],'classification':d['failure_classification_if_using_information_storage'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
