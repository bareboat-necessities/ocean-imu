"""Highest current finite prediction runtime entry for ALT.

Composes physical mean, attitude F/Q, integrated-OU Qaxis, BA decay and full
21-state covariance without accepting precomputed F_AA/Q_AA/F_LL/Q_LL/Q_BB or
phi factors. The shipping-level source entry additionally accepts a RawImuSample:
the wrapper's raw B-frame gyro is checked, then its exact shipping de-heel B->B'
map supplies the internal gyro used by nominal and covariance prediction.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as QX
from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as QB
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class QAxisBranch:
    correlated: bool
    sigma_aw: tuple
    marginal_psd: tuple
    final_psd: tuple
    machine_epsilon: F
    coefficient_branch: QB.Branch | None = None
    def __post_init__(self):
        if not isinstance(self.correlated,bool): raise TypeError('literal correlated/independent branch required')
        sig=M.mat(self.sigma_aw,3,3); eps=P.rational(self.machine_epsilon)
        if sig != M.transpose(sig) or eps <= 0: raise ValueError('symmetric Sigma_aw and positive epsilon required')
        expected=1 if self.correlated else 3
        if len(self.marginal_psd)!=expected or len(self.final_psd)!=expected: raise ValueError('Qaxis witness count must match shipping branch')
        if any(not isinstance(x,QX.PSDWitness) for x in self.marginal_psd+self.final_psd): raise TypeError('PSD runtime witnesses required')
        if self.coefficient_branch is not None and not isinstance(self.coefficient_branch,QB.Branch):
            raise TypeError('binary32 Qaxis coefficient branch required')
        object.__setattr__(self,'sigma_aw',tuple(map(tuple,sig))); object.__setattr__(self,'machine_epsilon',eps)

    def matrices(self,ou):
        small=None
        if self.coefficient_branch is not None:
            if self.coefficient_branch.tau_input != ou.tau or self.coefficient_branch.h_input != ou.h:
                raise ValueError('Qaxis binary32 formula branch detached from SAME OU tau/h operands')
            small=self.coefficient_branch.small
        if self.correlated:
            unit=QX.qaxis4(ou.tau,ou.h,1,ou.alpha,marginal_psd=self.marginal_psd[0],final_psd=self.final_psd[0],machine_epsilon=self.machine_epsilon,small_branch=small)
            return dict(qaxis_unit=unit,sigma_aw=self.sigma_aw,independent_qaxis=None)
        qs=[]
        for axis in range(3):
            qs.append(QX.qaxis4(ou.tau,ou.h,self.sigma_aw[axis][axis],ou.alpha,
                                marginal_psd=self.marginal_psd[axis],final_psd=self.final_psd[axis],machine_epsilon=self.machine_epsilon,small_branch=small))
        return dict(qaxis_unit=None,sigma_aw=None,independent_qaxis=qs)


def prediction(state,segment,*,gyro_body,angular,Qbase,ou,bias,qaxis:QAxisBranch,
               use_exact_attitude_Q=True,attitude_first_ldlt_success=True,
               attitude_second_ldlt_success=None):
    """Finite conditional predictor; gyro_body is the already de-heeled B' gyro."""
    if not isinstance(qaxis,QAxisBranch): raise TypeError('Qaxis runtime branch required')
    if angular.h != ou.h or ou.h != segment.h: raise ValueError('attitude, OU and physical segment must share one step')
    g=P.vec(gyro_body,3); ref=state.reference
    omega_hat=[g[i]-(ref.gyro_bias[i]-state.z[3+i]) for i in range(3)]
    if list(angular.w) != omega_hat: raise ValueError('attitude covariance angular rate detached from SAME nominal gyro prediction')
    qkw=qaxis.matrices(ou)
    return ATT.runtime_paired_prediction(
        state,segment,gyro_body=gyro_body,angular=angular,Qbase=Qbase,ou=ou,bias=bias,
        use_exact_Q=use_exact_attitude_Q,first_ldlt_success=attitude_first_ldlt_success,
        second_ldlt_success=attitude_second_ldlt_success,**qkw)


def prediction_from_raw(state,segment,sample:SENSOR.RawImuSample,*,angular,Qbase,ou,bias,qaxis:QAxisBranch,
                        use_exact_attitude_Q=True,attitude_first_ldlt_success=True,
                        attitude_second_ldlt_success=None):
    """Shipping-level prediction: raw B gyro -> same D_h -> internal B' predictor."""
    if not isinstance(sample,SENSOR.RawImuSample): raise TypeError('RawImuSample required')
    if sample.physical != segment.before:
        raise ValueError('raw IMU physical predecessor detached from prediction segment')
    g_internal=SENSOR.assert_prediction_gyro(sample,state,sample.raw_gyro_body)
    expected=sample.required_bias_corrected_relation(state.z[3:6])
    if list(angular.w) != list(expected):
        raise ValueError('attitude runtime omega detached from omega_sample+e_bg+n_g after de-heel')
    return prediction(state,segment,gyro_body=g_internal,angular=angular,Qbase=Qbase,ou=ou,bias=bias,qaxis=qaxis,
                      use_exact_attitude_Q=use_exact_attitude_Q,
                      attitude_first_ldlt_success=attitude_first_ldlt_success,
                      attitude_second_ldlt_success=attitude_second_ldlt_success)


def readiness():
    qb=QB.readiness()
    return {
      'free_prediction_transition_matrices_at_entry':False,
      'physical_mean_and_covariance_same_step':True,
      'same_bias_corrected_gyro_for_nominal_and_covariance':True,
      'raw_body_gyro_source_provenance_attached_at_shipping_entry':True,
      'shipping_deheel_map_before_internal_prediction_attached':True,
      'omega_sample_e_bg_n_g_relation_attached':True,
      'attitude_F_Q_runtime_graph_composed':True,
      'OU_mean_BA_runtime_graph_composed':True,
      'Qaxis_runtime_graph_composed':True,
      'Qaxis_binary32_formula_branch_graph_available':qb['small_general_comparison_binary32_attached'],
      'runtime_exp_trig_values_source_attached':False,
      'deheel_sincos_binary32_ancestry_attached':False,
      'sensor_residual_source_bounds_attached':False,
      'PSD_LDLT_eigensolver_finite_precision_attached':False,
      'same_history_tuner_parameters_attached':False,
      'post_prediction_prefixes_attached_separately':True,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
