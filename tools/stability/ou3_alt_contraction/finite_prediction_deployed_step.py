"""Bind the exact finite predictor to shipping's deployed small-step quaternion.

`prediction_quaternion_branch` proves that every nominal and shadow prediction
increment in the declared regional Live domain is strictly inside the 0.01 rad
polynomial branch.  Therefore no free quaternion coefficient remains here:
q_nominal and q_shadow are generated from the same omega_hat, e_bg and h that
feed the prediction.

This is exact real arithmetic. Binary32 FMA and normalization residuals remain
separate deployment supplies.
"""
from __future__ import annotations
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as G
from tools.stability.ou3_alt_contraction import prediction_quaternion_branch as BRANCH


def step_quaternion_polynomial(dtheta):
    d=G.vec(dtheta,3);u=sum((x*x for x in d),F(0))
    if u>=F(1,10000): raise ValueError('shipping polynomial branch requires ||dtheta|| < 0.01')
    u2=u*u
    w=F(1)-u/F(8)+u2/F(384)
    k=F(1,2)-u/F(48)+u2/F(3840)
    return [w,*[k*x for x in d]]


def prediction(mode,z,*,omega_hat,h,q15,coeff,w_bias,phi_true,phi_hat=None):
    z=G.vec(z,24);omega=G.vec(omega_hat,3);dt=G.rational(h)
    d_nom=[-omega[i]*dt for i in range(3)]
    d_shadow=[(-omega[i]+z[3+i])*dt for i in range(3)]
    qn=step_quaternion_polynomial(d_nom);qs=step_quaternion_polynomial(d_shadow)
    return G.prediction(mode,z,q_shadow=qs,q_nominal=qn,q15=q15,coeff=coeff,w_bias=w_bias,phi_true=phi_true,phi_hat=phi_hat)


def readiness():
    branch=BRANCH.build();fail=BRANCH.validate(branch)
    return {
      'qualification':'OU3_ALT_FINITE_DEPLOYED_PREDICTION_STEP_V1',
      'finite_prediction_graph_consumed':True,
      'source_uniform_polynomial_quaternion_branch_consumed':not fail,
      'same_omega_bg_h_generate_nominal_and_shadow_steps':True,
      'free_step_quaternion_coefficients_used':False,
      'trigonometric_branch_needed_in_regional_domain':False if not fail else None,
      'q15_same_history_source_relation_attached_to_finite_word':False,
      'bias_root_driver_temporal_relation_attached_to_finite_word':False,
      'covariance_frontend_successor_attached':False,
      'finite_precision_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
