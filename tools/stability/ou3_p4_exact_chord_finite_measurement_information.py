#!/usr/bin/env python3
"""Direct finite-state vector-chord + a_w information lower for P4.

This lemma is intentionally local to the selected complete-BRMM measurement
stack.  It upgrades the finite-angle *measurement* information from a
Jacobian/differential statement to an exact finite-residual statement.  It does
NOT by itself prove finite prediction/reset transport from a word-entry error.

Let p be the stacked tangent vector contribution from eta6 and let q be the
stacked exact finite-angle chord contribution.  The exact Cayley chord sector
holds component-eventwise and hence after stacking:

    q^T p = q^T q,                 ||q||^2 >= k ||p||^2.

Let u=C w be the selected accelerometer a_w contribution, with
||C||^2 <= c2, and let the same four actual S=0 records provide the directional
regularizer d_w ||w||^2.  The exact selected residual energy is

    ||q + u||^2 + d_w ||w||^2.

For any epsilon in (0,1), Young's inequality gives

    ||q+u||^2 >= (1-epsilon)||q||^2
                 -(1/epsilon-1)||u||^2.

Using ||p||^2 >= alpha ||eta6||^2, the chord lower and ||u||^2<=c2||w||^2,
we obtain a positive 2-block lower whenever epsilon>c2/(c2+d_w).  Eliminating
epsilon with the determinant/trace lower yields the simple rigorous bound

    lambda_finite >= (k*alpha*d_w)/(k*alpha + c2 + d_w).

This is the same algebraic form used by the finite-angle H18 backbone after
alpha -> k*alpha, but here it is proved directly for exact finite residuals;
there is no Taylor eta, no S^-1 scalarization and no packet-count factor.
Remaining non-a_w translation directions retain their directional four-S lower.

The result therefore closes the finite-state MEASUREMENT-STACK information
lemma.  Prediction/source transport, exact finite reset, BIAS1/projection and
prefix hard-domain retention remain downstream obligations and all P4 flags
stay false here.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_brmm_h18_information_composition as HINFO
import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD
import ou3_p4_finite_angle_h18_universal as H18FA

QUALIFICATION='OU3_P4_EXACT_CHORD_FINITE_MEASUREMENT_INFORMATION_V1'

def down(x):return math.nextafter(float(x),-math.inf)
def up(x):return math.nextafter(float(x),math.inf)

def build():
    h=HINFO.build();c=CHORD.build();fa=H18FA.build()
    bad={'hinfo':HINFO.validate(h),'chord':CHORD.validate(c),'finite_angle_h18':H18FA.validate(fa)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('finite chord measurement prerequisites failed: '+repr(bad))
    if not all(x['canonical_source']=='COMPLETE_BRMM_NORMAL_LIVE_WORD' for x in (h,c,fa)):
        raise RuntimeError('finite chord measurement lemma detached from complete BRMM')
    k=float(c['information_retention_factor_lower_full_entry'])
    alpha=float(h['eta6_information_lower'])
    comp=h['triangular_information_composition']
    d=float(comp['aw_direction_information_lower'])
    c2=float(comp['C_aw_spectral_norm_squared_upper'])
    a=down(k*alpha)
    trace=up(a+d+c2)
    lam=down(a*d/trace)
    nonaw=float(comp['non_aw_translation_lambda_min_lower'])
    D=down(min(lam,nonaw))
    if not all(math.isfinite(x) and x>0 for x in (k,alpha,d,c2,a,trace,lam,nonaw,D)):
        raise RuntimeError('finite chord measurement information lost positivity')
    # The differential backbone used precisely the same k*alpha replacement.
    # Agreement is a structural cross-check, not a reason to call the whole
    # finite map closed.
    fa_lam=float(fa['finite_angle_H18_information_lower'])
    rel=abs(D-fa_lam)/max(D,fa_lam)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'exact_finite_chord_sector_consumed':True,
      'finite_residual_not_Jacobian_only':True,
      'measurement_stack_only_scope':True,
      'selected_vector_PE_and_four_S_from_same_complete_word':True,
      'actual_applied_SpectralMSE_RS_retained_by_upstream_four_S':True,
      'exact_chord_information_factor_lower':k,
      'tangent_eta6_information_lower':alpha,
      'finite_eta6_information_lower':a,
      'aw_direction_information_lower':d,
      'accelerometer_aw_cross_norm_squared_upper':c2,
      'finite_coupled_eta6_aw_information_lower':lam,
      'non_aw_translation_information_lower':nonaw,
      'finite_H18_measurement_stack_information_lower':D,
      'differential_backbone_H18_information_lower_crosscheck':fa_lam,
      'finite_vs_differential_relative_difference':rel,
      'standalone_eta_Rinv_budget_used':False,
      'packet_count_multiplier_used':False,
      'independent_K_box_used':False,
      'prediction_transport_closed_here':False,
      'finite_reset_transport_closed_here':False,
      'source_uniform_endpoint_augmented_LDLT_closed_here':False,
      'source_uniform_every_prefix_augmented_LDLT_closed_here':False,
      'same_graph_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for x in ('exact_finite_chord_sector_consumed','finite_residual_not_Jacobian_only','measurement_stack_only_scope',
              'selected_vector_PE_and_four_S_from_same_complete_word','actual_applied_SpectralMSE_RS_retained_by_upstream_four_S'):
        if d.get(x) is not True:f.append(x+' not true')
    for x in ('standalone_eta_Rinv_budget_used','packet_count_multiplier_used','independent_K_box_used',
              'prediction_transport_closed_here','finite_reset_transport_closed_here','source_uniform_endpoint_augmented_LDLT_closed_here',
              'source_uniform_every_prefix_augmented_LDLT_closed_here','same_graph_every_prefix_hard_domain_retention_closed_here',
              'P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(x) is not False:f.append(x+' not false')
    for x in ('exact_chord_information_factor_lower','finite_eta6_information_lower','finite_coupled_eta6_aw_information_lower',
              'finite_H18_measurement_stack_information_lower'):
        if not float(d.get(x,0))>0:f.append(x+' not positive')
    if float(d.get('finite_vs_differential_relative_difference',math.inf))>5e-12:
        f.append('finite direct measurement lower disagrees with differential algebra crosscheck')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'k':d['exact_chord_information_factor_lower'],'finite_eta6':d['finite_eta6_information_lower'],
      'finite_coupled':d['finite_coupled_eta6_aw_information_lower'],'D18':d['finite_H18_measurement_stack_information_lower'],
      'crosscheck_rel':d['finite_vs_differential_relative_difference'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
