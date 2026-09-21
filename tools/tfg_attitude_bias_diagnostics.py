#!/usr/bin/env python3
"""Static structural audit of TFG versus OU-III accelerometer attitude/bias update.

This study intentionally avoids replaying the long tuning campaigns. It records
the exact physical residual/Jacobian conventions and the covariance/reset
operations that can redistribute an accelerometer innovation among tilt, a_w
and b_a. Numerical replay diagnostics are produced separately.
"""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parents[1]
TFG=ROOT/'src/kalman_tfg/Kalman3D_Wave_TFG.h'
OU3=ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
def grab(text,start,end):
 a=text.index(start);b=text.index(end,a);return text[a:b]
def main():
 t=TFG.read_text();o=OU3.read_text()
 ts=grab(t,'    void accel_residual','    void mag_residual')
 tu=grab(t,'    bool apply_update3_','    void reorthonormalize_')
 os=grab(o,'void Kalman3D_Wave_OU_III<T, with_gyro_bias, with_accel_bias>::measurement_update_acc_only','template<typename T, bool with_gyro_bias, bool with_accel_bias>\nvoid Kalman3D_Wave_OU_III<T, with_gyro_bias, with_accel_bias>::measurement_update_mag_only')
 checks={
  'tfg_world_residual': 'X_.R * (acc_meas_body - ba)' in ts,
  'tfg_H_theta_gravity_only': 'skew<T>(gravity_world())' in ts,
  'tfg_H_aw_minus_I': 'OFF_AW) = -Matrix3::Identity()' in ts,
  'tfg_H_ba_level2_minus_I': 'Matrix3(-Matrix3::Identity())' in ts,
  'tfg_full_joseph': 'IKH * P_ * IKH.transpose() + K * Rw * K.transpose()' in tu,
  'tfg_full_group_reset': 'build_reset_jacobian(correction, Jreset)' in tu,
  'ou3_body_residual': 'f_meas - f_pred' in os,
  'ou3_H_theta_specific_force': 'J_att = -skew_symmetric_matrix(f_cog_b)' in os,
  'ou3_H_aw_Rwb': 'J_aw  =  R_wb()' in os,
  'ou3_H_ba_plus_I': 'PCt.noalias() += P_all_ba' in os,
  'ou3_quaternion_reset': 'applyQuaternionCorrectionFromErrorState()' in os,
 }
 if not all(checks.values()):raise SystemExit('structural signature changed: '+repr(checks))
 finding={
  'checks':checks,
  'principal_difference':[
   'TFG forms the accelerometer innovation in WORLD coordinates after rotating the bias-corrected body measurement.',
   'Its invariant attitude Jacobian is skew(g), while a_w and level-2 b_a enter as -I.',
   'OU-III forms the innovation in BODY coordinates; its attitude Jacobian is -skew(R_wb(a_w-g)), a_w enters through R_wb, and b_a through +I.',
   'TFG applies the full TwoFrameGroup finite-retraction covariance reset after every accepted accelerometer update; OU-III applies its MEKF quaternion correction/reset.',
  ],
  'interpretation':'These are coordinate-consistent formulations in principle, so different Jacobian appearance alone is not an error. The replay study must test physical innovation allocation and covariance transport before/after reset.',
 }
 out=ROOT/'reports/results/tfg_attitude_bias_diagnostics';out.mkdir(parents=True,exist_ok=True)
 (out/'structural-audit.json').write_text(json.dumps(finding,indent=2)+'\n')
 print(json.dumps(finding,indent=2))
if __name__=='__main__':main()
