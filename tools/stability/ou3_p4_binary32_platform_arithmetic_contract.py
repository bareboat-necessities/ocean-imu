#!/usr/bin/env python3
"""Conditional binary32 execution contract for non-source library arithmetic.

The explicit OU-III Joseph/reset scalar loops are enclosed independently. ISO
C++ does not guarantee Eigen LDLT or libm accuracy, so the remaining numerical
facts are explicit execution premises, parallel to canonical P3 execution
premises. Deployment qualification remains separate.

Two directly consumable forward postconditions are supplied:

* gain rows after safe_ldlt3_/Eigen solve; and
* the complete homogeneous time-update state map after all libm-generated OU
  and rotation coefficients.

The second postcondition is preferable to propagating a raw libm ULP bound
through cancellation-sensitive coefficient formulas. The target qualifier may
use the ULP contract as implementation evidence, but P4 consumes the complete
map postcondition.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

U=2.0**-24
LDLT_RESIDUAL_REL_UPPER=math.nextafter(256.0*U/(1.0-256.0*U),math.inf)
LDLT_GAIN_ROW_FORWARD_REL_UPPER=2.0**-12
LIBM_ULP_UPPER=8.0
# ||F_fp(x)-F_real(x)|| <= eps_time*max(1,||x||), over one admitted IMU prediction.
TIME_UPDATE_STATE_FORWARD_REL_UPPER=2.0**-11


def build()->dict:
    return {
      'qualification':'OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V3',
      'runtime_scalar':'IEEE754_binary32','rounding_mode_required':'roundTiesToEven',
      'finite_inputs_required':True,'flush_to_zero_must_not_affect_admitted_nonzero_arguments':True,
      'compiler_fast_math_must_not_break_declared_postconditions':True,'unit_roundoff':U,
      'eigen_ldlt3_solve_residual_postcondition':'||S*xhat-r||_2 <= eps_solve*(||S||_2*||xhat||_2+||r||_2)',
      'eigen_ldlt3_relative_residual_upper':LDLT_RESIDUAL_REL_UPPER,
      'eigen_ldlt3_gain_row_forward_postcondition':'||Khat_i-K_i||_2 <= eps_gain*max(1,||K_i||_2)',
      'eigen_ldlt3_gain_row_forward_relative_upper':LDLT_GAIN_ROW_FORWARD_REL_UPPER,
      'gain_forward_postcondition_is_platform_premise_not_condition_number_inference':True,
      'safe_ldlt3_diagonal_bump_branch_included':True,
      'libm_functions':['sin','cos','exp','expm1','sqrt'],
      'libm_finite_result_ulp_error_upper':LIBM_ULP_UPPER,
      'libm_argument_domain_is_existing_Normal_Live_compact_domain':True,
      'time_update_state_forward_postcondition':'||F_fp(x)-F_real(x)||_2 <= eps_time*max(1,||x||_2) per admitted IMU prediction',
      'time_update_state_forward_relative_upper':TIME_UPDATE_STATE_FORWARD_REL_UPPER,
      'time_update_forward_postcondition_includes_libm_and_derived_coefficient_arithmetic':True,
      'source_domain_changed':False,'filter_changed':False,'quality_gate_changed':False,
      'conditional_mathematical_execution_premise':True,
      'target_toolchain_qualified_here':False,'assembled_device_qualified_here':False,
      'deployment_claim_authorized_here':False,
      'P4_may_consume_as_explicit_execution_premise':True,'P4_DEPLOYMENT_PASS':False,
    }


def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V3':f.append('qualification mismatch')
    for k in ('finite_inputs_required','flush_to_zero_must_not_affect_admitted_nonzero_arguments',
              'compiler_fast_math_must_not_break_declared_postconditions','gain_forward_postcondition_is_platform_premise_not_condition_number_inference',
              'safe_ldlt3_diagonal_bump_branch_included','libm_argument_domain_is_existing_Normal_Live_compact_domain',
              'time_update_forward_postcondition_includes_libm_and_derived_coefficient_arithmetic',
              'conditional_mathematical_execution_premise','P4_may_consume_as_explicit_execution_premise'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_domain_changed','filter_changed','quality_gate_changed','target_toolchain_qualified_here',
              'assembled_device_qualified_here','deployment_claim_authorized_here','P4_DEPLOYMENT_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    er=float(d.get('eigen_ldlt3_relative_residual_upper',math.nan));ef=float(d.get('eigen_ldlt3_gain_row_forward_relative_upper',math.nan));et=float(d.get('time_update_state_forward_relative_upper',math.nan));ulp=float(d.get('libm_finite_result_ulp_error_upper',math.nan))
    if not (math.isfinite(er) and 0<er<1e-3):f.append('invalid LDLT residual premise')
    if not (math.isfinite(ef) and 0<ef<1e-2):f.append('invalid LDLT forward premise')
    if not (math.isfinite(et) and 0<et<1e-2):f.append('invalid time-update forward premise')
    if not (math.isfinite(ulp) and 1<=ulp<=16):f.append('invalid libm ULP premise')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'ldlt_residual':d['eigen_ldlt3_relative_residual_upper'],'ldlt_forward':d['eigen_ldlt3_gain_row_forward_relative_upper'],'time_forward':d['time_update_state_forward_relative_upper'],'libm_ulp':d['libm_finite_result_ulp_error_upper'],'deployment':d['P4_DEPLOYMENT_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
