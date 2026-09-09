#!/usr/bin/env python3
"""Conditional binary32 execution contract for the two non-source arithmetic pieces.

P4 already encloses explicit scalar loops in the OU-III Joseph/reset path.  Two
implementation-dependent operations cannot be inferred from ISO C++ alone:

* Eigen's fixed-size 3x3 LDLT factorization/solve used for each accepted vector
  correction; and
* the target libm implementations of sin/cos/exp/expm1/sqrt used to construct
  bounded time-update/reset coefficients.

This file makes those requirements *execution premises*, exactly like the
canonical P3 execution premises.  It is deliberately not a hardware/device
qualification.  The mathematical P4 implication may consume the premise; a
shipping/deployment claim additionally has to qualify the target toolchain.

The contract is stated as residual/ULP postconditions, not as a claim about a
particular Eigen or libm internal algorithm:

  ||S xhat-r||_2 <= eps_solve (||S||_2 ||xhat||_2 + ||r||_2)

for each 3x3 solve after the shipping safe_ldlt3_ branch, and a uniform ULP
error for finite libm outputs on the already-declared compact Normal-Live
argument ranges.  This permits the downstream finite-precision producer to
turn both effects into an additive ISS input without changing the estimator.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

U = 2.0**-24
# Deliberately generous postcondition constants.  They are premises to qualify,
# not empirical measurements and not silently attributed to the C++ standard.
LDLT_RESIDUAL_REL_UPPER = math.nextafter(256.0 * U / (1.0 - 256.0 * U), math.inf)
LIBM_ULP_UPPER = 8.0


def build() -> dict:
    return {
        "qualification": "OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V1",
        "runtime_scalar": "IEEE754_binary32",
        "rounding_mode_required": "roundTiesToEven",
        "finite_inputs_required": True,
        "flush_to_zero_must_not_affect_admitted_nonzero_arguments": True,
        "compiler_fast_math_must_not_break_declared_postconditions": True,
        "unit_roundoff": U,
        "eigen_ldlt3_solve_postcondition": (
            "||S*xhat-r||_2 <= eps_solve*(||S||_2*||xhat||_2+||r||_2)"
        ),
        "eigen_ldlt3_relative_residual_upper": LDLT_RESIDUAL_REL_UPPER,
        "safe_ldlt3_diagonal_bump_branch_included": True,
        "libm_functions": ["sin", "cos", "exp", "expm1", "sqrt"],
        "libm_finite_result_ulp_error_upper": LIBM_ULP_UPPER,
        "libm_argument_domain_is_existing_Normal_Live_compact_domain": True,
        "source_domain_changed": False,
        "filter_changed": False,
        "quality_gate_changed": False,
        "conditional_mathematical_execution_premise": True,
        "target_toolchain_qualified_here": False,
        "assembled_device_qualified_here": False,
        "deployment_claim_authorized_here": False,
        "P4_may_consume_as_explicit_execution_premise": True,
        "P4_DEPLOYMENT_PASS": False,
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("qualification") != "OU3_P4_CONDITIONAL_BINARY32_PLATFORM_ARITHMETIC_V1":
        f.append("qualification mismatch")
    for k in (
        "finite_inputs_required", "flush_to_zero_must_not_affect_admitted_nonzero_arguments",
        "compiler_fast_math_must_not_break_declared_postconditions",
        "safe_ldlt3_diagonal_bump_branch_included",
        "libm_argument_domain_is_existing_Normal_Live_compact_domain",
        "conditional_mathematical_execution_premise",
        "P4_may_consume_as_explicit_execution_premise",
    ):
        if d.get(k) is not True: f.append(k+" not true")
    for k in (
        "source_domain_changed", "filter_changed", "quality_gate_changed",
        "target_toolchain_qualified_here", "assembled_device_qualified_here",
        "deployment_claim_authorized_here", "P4_DEPLOYMENT_PASS",
    ):
        if d.get(k) is not False: f.append(k+" not false")
    e=float(d.get("eigen_ldlt3_relative_residual_upper", math.nan))
    if not (math.isfinite(e) and 0.0 < e < 1e-3): f.append("invalid LDLT residual premise")
    ulp=float(d.get("libm_finite_result_ulp_error_upper", math.nan))
    if not (math.isfinite(ulp) and 1.0 <= ulp <= 16.0): f.append("invalid libm ULP premise")
    return f


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"ldlt_eps":d["eigen_ldlt3_relative_residual_upper"],"libm_ulp":d["libm_finite_result_ulp_error_upper"],"deployment":d["P4_DEPLOYMENT_PASS"],"failures":f},sort_keys=True))
    return int(bool(f))

if __name__ == "__main__": raise SystemExit(main())
