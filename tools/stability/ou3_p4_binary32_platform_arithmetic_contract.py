#!/usr/bin/env python3
"""Conditional binary32 execution contract for non-source library arithmetic.

The explicit OU-III Joseph/reset scalar loops are enclosed independently.  ISO
C++ alone does not provide accuracy contracts for Eigen's fixed-size LDLT solve
or the target libm.  P4 therefore states those two facts as execution premises,
parallel to the existing conditional P3 execution premises; target deployment
qualification remains separate.

The LDLT contract includes both a residual check and a directly consumable
forward gain-row postcondition.  The latter is intentional: a gross
source-uniform condition-number ceiling is too pessimistic to turn the residual
bound into a useful forward bound.  A target qualification must check the
forward postcondition on the declared compact innovation family; the proof does
not infer it from a replay or from ISO C++.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

U=2.0**-24
LDLT_RESIDUAL_REL_UPPER=math.nextafter(256.0*U/(1.0-256.0*U),math.inf)
# Relative row error against max(1,||K_i||).  This is an execution premise,
# deliberately loose compared with ordinary binary32 behavior.
LDLT_GAIN_ROW_FORWARD_REL_UPPER=2.0**-12
LIBM_ULP_UPPER=8.0


def build()->dict:
    return {
      'qualification':'OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V2',
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
      'source_domain_changed':False,'filter_changed':False,'quality_gate_changed':False,
      'conditional_mathematical_execution_premise':True,
      'target_toolchain_qualified_here':False,'assembled_device_qualified_here':False,
      'deployment_claim_authorized_here':False,
      'P4_may_consume_as_explicit_execution_premise':True,'P4_DEPLOYMENT_PASS':False,
    }


def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V2':f.append('qualification mismatch')
    for k in ('finite_inputs_required','flush_to_zero_must_not_affect_admitted_nonzero_arguments',
              'compiler_fast_math_must_not_break_declared_postconditions','gain_forward_postcondition_is_platform_premise_not_condition_number_inference',
              'safe_ldlt3_diagonal_bump_branch_included','libm_argument_domain_is_existing_Normal_Live_compact_domain',
              'conditional_mathematical_execution_premise','P4_may_consume_as_explicit_execution_premise'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_domain_changed','filter_changed','quality_gate_changed','target_toolchain_qualified_here',
              'assembled_device_qualified_here','deployment_claim_authorized_here','P4_DEPLOYMENT_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    er=float(d.get('eigen_ldlt3_relative_residual_upper',math.nan)); ef=float(d.get('eigen_ldlt3_gain_row_forward_relative_upper',math.nan)); ulp=float(d.get('libm_finite_result_ulp_error_upper',math.nan))
    if not (math.isfinite(er) and 0<er<1e-3):f.append('invalid LDLT residual premise')
    if not (math.isfinite(ef) and 0<ef<1e-2):f.append('invalid LDLT forward premise')
    if not (math.isfinite(ulp) and 1<=ulp<=16):f.append('invalid libm ULP premise')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'ldlt_residual':d['eigen_ldlt3_relative_residual_upper'],'ldlt_forward':d['eigen_ldlt3_gain_row_forward_relative_upper'],'libm_ulp':d['libm_finite_result_ulp_error_upper'],'deployment':d['P4_DEPLOYMENT_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
