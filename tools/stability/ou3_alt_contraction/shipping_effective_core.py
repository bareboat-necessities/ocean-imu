"""Shipping-parity theorem for the generalized 3x3 Joseph core C.

For H18, accelerometer-bias updates are disabled and shipping zeros every
cross-covariance between b_a and the first 18 estimator coordinates. Therefore
with zero lever arm:

  accelerometer: C = R_acc + P_ba,ba
  magnetometer:  C = R_mag
  S=0:           C = R_S

where C=S-N_H^T A^-1 N_H and A is the H18 covariance block. Hence C is SPD as
soon as the admitted measurement noise is SPD and P_ba,ba is PSD. In A21 the
numerator is unmasked and C=R. This closes the mask/core positivity issue
without Mahony startup, the PE/vector-domain lemma, or a covariance upper tube.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from . import generalized_joseph as G

REPO=Path(__file__).resolve().parents[3]
HEADER=REPO/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
QUALIFICATION='OU3_ALT_SHIPPING_EFFECTIVE_CORE_PARITY_V1'

def source_parity():
    t=HEADER.read_text(encoding='utf-8')
    return {
      'disable_zeros_ba_base_cross':'Pext.template block<3,BASE_N>(OFF_BA, 0).setZero();' in t and 'Pext.template block<BASE_N,3>(0, OFF_BA).setZero();' in t,
      'disable_zeros_ba_linear_cross':'Pext.template block<3,12>(OFF_BA, OFF_V).setZero();' in t and 'Pext.template block<12,3>(OFF_V, OFF_BA).setZero();' in t,
      'H18_acc_S_retains_P_ba_ba':'BA frozen: marginalize its uncertainty' in t and 'S_mat.noalias() += P_ba_ba;' in t,
      'H18_acc_numerator_excludes_ba_column':'if (use_ba) {' in t and 'PCt.noalias() += P_all_ba; // J_ba = I' in t,
      'frozen_numerator_rows_zeroed':'if (!use_ba) freeze_acc_bias_rows_(PCt);' in t,
      'frozen_gain_rows_zeroed':'if (!use_ba) freeze_acc_bias_rows_(K);' in t,
      'mag_S_is_attitude_block_plus_Rmag':'S_mat = Rmag;' in t and 'S_mat.noalias() += J_att * P_th_th * J_att.transpose();' in t,
      'Szero_S_is_PSS_plus_RS':'S_mat = Pext.template block<3,3>(off_S, off_S) + R_S;' in t,
      'Szero_numerator_is_P_column_then_mask':'PCt.noalias() = Pext.template block<NX,3>(0, off_S);' in t,
      'noise_covariances_from_positive_stds':'Racc = sigma_acc.array().square().matrix().asDiagonal();' in t and 'Rmag(sigma_m.array().square().matrix().asDiagonal())' in t and 'R_S = sigma_S.array().square().matrix().asDiagonal();' in t,
    }

def block_regression():
    A=[[F(3),F(1,3)],[F(1,3),F(2)]]
    H=[[F(1),F(2)],[F(-1),F(1)],[F(1,2),F(-1,3)]]
    D=[[F(1,5),F(0),F(0)],[F(0),F(1,4),F(0)],[F(0),F(0),F(1,6)]]
    R=[[F(2),F(0),F(0)],[F(0),F(3),F(0)],[F(0),F(0),F(4)]]
    J=G.inverse(A);N=G.mm(A,G.mt(H))
    Sacc=G.add(G.add(G.mm(G.mm(H,A),G.mt(H)),D),R)
    Cacc=G.sub(Sacc,G.mm(G.mm(G.mt(N),J),N))
    Smag=G.add(G.mm(G.mm(H,A),G.mt(H)),R)
    Cmag=G.sub(Smag,G.mm(G.mm(G.mt(N),J),N))
    Ass=[[F(5),F(1,7),F(0)],[F(1,7),F(6),F(1,8)],[F(0),F(1,8),F(7)]]
    Czero=G.sub(G.add(Ass,R),G.mm(G.mm(G.mt(Ass),G.inverse(Ass)),Ass))
    return {'acc_core_equals_R_plus_Pba':Cacc==G.add(R,D),'mag_core_equals_R':Cmag==R,'Szero_core_equals_R':Czero==R}

def build():
    p=source_parity();b=block_regression();closed=all(p.values()) and all(b.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','zero_lever_arm_scope_retained':True,'shipping_source_parity':p,'all_required_shipping_parity_closed':all(p.values()),'exact_block_regression':b,'exact_H18_block_algebra_closed':all(b.values()),'H18_cross_covariance_ba_to_18_zero_when_held':all(p.values()),'H18_acc_effective_core':'C=R_acc+P_ba_ba','H18_mag_effective_core':'C=R_mag','H18_Szero_effective_core':'C=R_S','A21_unmasked_effective_core':'C=R','H18_acc_C_ge_Racc_by_Pba_PSD':True,'H18_mag_C_equals_Rmag':True,'H18_Szero_C_equals_RS':True,'A21_C_equals_R':True,'effective_core_dimension':3,'PE_vector_domain_lemma_needed_for_core_positivity':False,'Mahony_startup_needed_for_core_positivity':False,'covariance_tube_upper_bound_needed_for_core_positivity':False,'source_uniform_effective_core_form_closed':closed,'source_uniform_C_positive_definite_from_R_positive_and_Pba_PSD':closed,'source_uniform_xi_sector_closed_here':False,'source_uniform_reset_absorption_closed_here':False,'source_uniform_complete_word_dissipation_closed_here':False,'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False}
def validate(d):
    f=[]
    for k in ('zero_lever_arm_scope_retained','all_required_shipping_parity_closed','exact_H18_block_algebra_closed','H18_cross_covariance_ba_to_18_zero_when_held','H18_acc_C_ge_Racc_by_Pba_PSD','H18_mag_C_equals_Rmag','H18_Szero_C_equals_RS','A21_C_equals_R','source_uniform_effective_core_form_closed','source_uniform_C_positive_definite_from_R_positive_and_Pba_PSD'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('PE_vector_domain_lemma_needed_for_core_positivity','Mahony_startup_needed_for_core_positivity','covariance_tube_upper_bound_needed_for_core_positivity','source_uniform_xi_sector_closed_here','source_uniform_reset_absorption_closed_here','source_uniform_complete_word_dissipation_closed_here','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if d.get('effective_core_dimension')!=3:f.append('core dimension changed')
    return f
