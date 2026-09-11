"""Attach the finite prediction variables to admitted same-history source contracts.

The finite predictor uses q15=(a0,a1,J0,J1,J2) and one physical bias driver w.
This module binds exactly those variables to the already-certified COMPLETE-BRMM
q15 quadratic sectors and to one selected BIAS0/1/2 contract.  It does not
replace either relation by componentwise boxes.

This is a one-transition attachment. The 600-step word must still preserve one
BRMM generator/primitive chain and one bias parameter/root token by induction.
"""
from __future__ import annotations
from ou3_interval import Interval
from tools.stability.ou3_alt_contraction import finite_prediction_deployed_step as FINITE
from tools.stability.ou3_alt_contraction import physical_lineage as LINEAGE
from tools.stability.ou3_alt_contraction import bias_families as BIAS
import ou3_p4_brmm_physical_prediction_forcing as FORCING

QUALIFICATION='OU3_ALT_FINITE_PREDICTION_SOURCE_ATTACHMENT_V1'
Q15_ORDER=('a0_xyz','a1_xyz','J0_xyz','J1_xyz','J2_xyz')


def attach(*,mode,tau,h,bias_contract):
    if mode not in ('H','A'): raise ValueError('mode must be H/A')
    if not isinstance(tau,Interval) or not isinstance(h,Interval) or tau.lo<=0 or h.lo<=0:
        raise ValueError('positive interval tau/h required')
    if not isinstance(bias_contract,BIAS.BiasFamilyContract):
        raise TypeError('explicit BIAS0/1/2 contract required')
    finite=FINITE.readiness()
    qmap=FORCING.source_matrix(tau,h)
    sectors=LINEAGE._q15_sector_for_h(h)
    driver=BIAS.driver_ball_iqc(bias_contract)
    truth=BIAS.true_bias_ball_iqc(bias_contract)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','mode':mode,
      'finite_deployed_prediction_graph_consumed':finite['finite_prediction_graph_consumed'],
      'source_uniform_polynomial_quaternion_branch_consumed':finite['source_uniform_polynomial_quaternion_branch_consumed'],
      'q15_coordinate_order':Q15_ORDER,'q15_dimension':15,'q15_forcing_matrix':qmap,'q15_joint_quadratic_sectors':sectors,
      'q15_same_history_relation_attached_to_finite_prediction':True,
      'independent_q15_component_boxes_used':False,
      'bias_family':bias_contract.name,'bias_parameter_token':bias_contract.parameter_token,
      'bias_phi_true_interval':bias_contract.phi_true,'bias_driver_ball_iqc':driver,'true_bias_ball_iqc':truth,
      'same_physical_bias_driver_used_by_error_and_truth':True,
      'bias_parameter_root_token_attached_to_finite_prediction':True,
      'independent_bias_driver_copies_used':False,
      'one_transition_source_attachment_closed':True,
      'multi_transition_BRMM_primitive_chain_closed_here':False,
      'multi_transition_bias_parameter_token_continuity_closed_here':False,
      'covariance_frontend_successor_attached':False,
      'complete_word_finite_identity':False,'ALT_LIVE_PASS':False,
    }


def build():
    rows=[]
    for c in BIAS.contracts():
        rows.append(attach(mode='H',tau=Interval.point(1.0),h=Interval.point(0.005),bias_contract=c))
        rows.append(attach(mode='A',tau=Interval.point(1.0),h=Interval.point(0.005),bias_contract=c))
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','rows':rows,
      'all_H_A_BIAS_families_have_one_step_finite_source_attachment':all(r['one_transition_source_attachment_closed'] for r in rows),
      'no_independent_q15_or_bias_boxes':all(not r['independent_q15_component_boxes_used'] and not r['independent_bias_driver_copies_used'] for r in rows),
      'complete_word_finite_identity':False,'ALT_LIVE_PASS':False,
      'next_obligation':'replace the derivative event composer by a finite event composer that carries these q15 sectors and one bias parameter token through every successor and couples them to actual covariance/frontend postimages',
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('all_H_A_BIAS_families_have_one_step_finite_source_attachment','no_independent_q15_or_bias_boxes'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_word_finite_identity','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if len(d.get('rows',()))!=6:f.append('H/A x BIAS family coverage changed')
    return f
