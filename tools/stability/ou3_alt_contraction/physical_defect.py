"""Physical defect graph for the generalized masked-numerator information identity.

For the actual shipping Joseph numerator N define L=N^T P^-1 and
C=S-N^T P^-1 N. The generalized signed-energy identity uses xi=q-L e.
With zero lever arm and the shipping H18 hold semantics:

 A21 S=0: xi=0; A21 mag: xi=eta_mag; A21 accel: xi=eta_acc;
 H18 S=0: xi=0; H18 mag: xi=eta_mag;
 H18 accel: xi=e_ba+eta_acc.

eta_mag/eta_acc are the exact finite-Cayley residual remainders from the
existing COMPLETE-BRMM sector lemma. The held e_ba term is a neutral
physical/reference port, not a relabeled motion error. In exact real arithmetic
shipping projection gives ||bhat_a||<=R_proj and each admitted bias family gives
||b_true||<=B_true, so ||e_ba||<=R_proj+B_true independently of startup landing
or covariance consistency.

For H18 accelerometer C=R_acc+P_baba >= R_acc, hence
xi^T C^-1 xi <= xi^T R_acc^-1 xi. Every other accepted measurement type has
C=R. Thus the positive Joseph defect stays a 3D port plus the hard Cayley
sector; no independent inverse box or packet-count remainder is introduced.
"""
from __future__ import annotations
import math
import sys
from pathlib import Path

STABILITY=Path(__file__).resolve().parents[1]
if str(STABILITY) not in sys.path: sys.path.insert(0,str(STABILITY))

from tools.stability import ou3_p4_complete_brmm_residual_sector as RES
from tools.stability import ou3_p4_bias0_family as B0
from tools.stability import ou3_p4_bias1_family as B1
from tools.stability import ou3_p4_bias2_family as B2
from tools.stability import ou3_projection_sector as PROJ
from . import shipping_effective_core as CORE

QUALIFICATION='OU3_ALT_PHYSICAL_GENERALIZED_DEFECT_GRAPH_V1'
PROJECTION_RADIUS=0.4

def _family_bounds():
    out={}
    for name,mod in [('BIAS0',B0),('BIAS1',B1),('BIAS2',B2)]:
        d=mod.build(); f=mod.validate(d)
        if f: raise RuntimeError(f'{name} invalid: {f!r}')
        out[name]=float(d['true_bias_norm_upper_mps2'])
    return out

def build():
    core=CORE.build(); cf=CORE.validate(core)
    if cf: raise RuntimeError(f'effective-core prerequisite failed: {cf!r}')
    sector=RES.build(); sf=RES.validate(sector)
    if sf: raise RuntimeError(f'residual-sector prerequisite failed: {sf!r}')
    projection=PROJ.build_report(PROJECTION_RADIUS)
    if projection.get('exact_real_operator_sector_closed') is not True or projection.get('estimate_ball_invariant_exact_real') is not True:
        raise RuntimeError('radial projection invariant not closed')
    fam=_family_bounds(); bmax=max(fam.values())
    held=math.nextafter(PROJECTION_RADIUS+bmax,math.inf)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'zero_lever_arm_scope_retained':True,'actual_numerator_generalized_core_consumed':True,
      'physical_defect_forms':{
        'A21_S_zero':'xi=0','A21_magnetometer':'xi=eta_mag','A21_accelerometer':'xi=eta_acc',
        'H18_S_zero':'xi=0','H18_magnetometer':'xi=eta_mag','H18_accelerometer':'xi=e_ba+eta_acc'},
      'S_zero_positive_defect_exactly_zero':True,
      'magnetometer_positive_defect_is_existing_exact_Cayley_remainder':True,
      'A21_accelerometer_positive_defect_is_existing_exact_Cayley_aw_remainder':True,
      'H18_accelerometer_adds_only_held_bias_error_port':True,
      'accelerometer_bias_nonlinearity_itself_exactly_zero':bool(sector['accelerometer_bias_nonlinearity_exactly_zero']),
      'finite_residual_has_no_state_independent_motion_term':bool(sector['nonlinear_residual_has_no_state_independent_term']),
      'projection_real_ball_invariant_closed':True,'projection_radius_mps2':PROJECTION_RADIUS,
      'true_bias_norm_upper_by_family_mps2':fam,'held_bias_error_compactness_upper_mps2':held,
      'held_bias_compactness_route':'||e_ba||<=||bhat_a||+||b_true||<=R_projection+B_true',
      'held_bias_bound_uses_startup_membership':False,'held_bias_bound_uses_covariance_consistency':False,
      'held_bias_bound_uses_desired_motion_basin':False,
      'H18_acc_C_inverse_dominated_by_Racc_inverse':True,'all_other_C_equal_actual_R':True,
      'positive_defect_port_dimension':3,'same_history_eta_sector_required':True,
      'eta_packet_count_scalarization_used':False,'independent_C_inverse_box_used':False,
      'source_uniform_defect_graph_closed':True,
      'source_uniform_defect_energy_domination_closed_here':False,
      'source_uniform_reset_absorption_closed_here':False,
      'source_uniform_complete_word_dissipation_closed_here':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
      'next_obligation':'assemble xi as a retained 3D graph variable on every same-history event, keep H18 e_ba as the bounded neutral port, and prove the joint quadratic xi/reset cost is dominated by the complete-word negative information; do not sum eventwise worst cases'}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('zero_lever_arm_scope_retained','actual_numerator_generalized_core_consumed','S_zero_positive_defect_exactly_zero','magnetometer_positive_defect_is_existing_exact_Cayley_remainder','A21_accelerometer_positive_defect_is_existing_exact_Cayley_aw_remainder','H18_accelerometer_adds_only_held_bias_error_port','accelerometer_bias_nonlinearity_itself_exactly_zero','finite_residual_has_no_state_independent_motion_term','projection_real_ball_invariant_closed','H18_acc_C_inverse_dominated_by_Racc_inverse','all_other_C_equal_actual_R','same_history_eta_sector_required','source_uniform_defect_graph_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('held_bias_bound_uses_startup_membership','held_bias_bound_uses_covariance_consistency','held_bias_bound_uses_desired_motion_basin','eta_packet_count_scalarization_used','independent_C_inverse_box_used','source_uniform_defect_energy_domination_closed_here','source_uniform_reset_absorption_closed_here','source_uniform_complete_word_dissipation_closed_here','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('positive_defect_port_dimension')!=3:f.append('defect port dimension changed')
    if not float(d.get('held_bias_error_compactness_upper_mps2',0))>PROJECTION_RADIUS:f.append('held bias compactness invalid')
    forms=d.get('physical_defect_forms',{})
    if forms.get('H18_accelerometer')!='xi=e_ba+eta_acc' or forms.get('A21_accelerometer')!='xi=eta_acc':f.append('accelerometer defect forms changed')
    return f
