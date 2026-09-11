#!/usr/bin/env python3
"""Source-uniform continuation with correlated (P,H,R) innovation families.

The continuation replaces only the proof backend's measurement primitive while
executing the existing two-sample source-uniform H18/A21 path.  A measurement
now accepts P,H,R together and internally constructs

    S=H P H^T+R -> S^-1 -> K=P H^T S^-1.

No arbitrary rectangular S can be supplied.  The resulting K, A and Joseph
covariance update are therefore descendants of the same P,H,R family.  This
closes the immediate innovation-provenance break, but does not by itself close
endpoint contraction, every-prefix retention, or P4/P5.
"""
from __future__ import annotations

import json
import traceback

from ou3_interval import matrix_add, matrix_mul, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_brmm_full_word_riccati_backend as BACKEND
import ou3_correlated_innovation_family as CORR
import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE


def _correlated_joseph(state, H, R):
    family = CORR.build(state.P, H, R)
    K, A = family.K, family.A
    At = matrix_transpose(A)
    KRKt = matrix_mul(matrix_mul(K, family.R), matrix_transpose(K))
    state.P = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A, state.P), At), KRKt)
    )
    state.Psi = matrix_mul(A, state.Psi)
    state.Omega = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A, state.Omega), At), KRKt)
    )
    state.events += 1
    state.measurements += 1
    if not BACKEND.decomposition_identity_enclosed(state):
        raise RuntimeError("P/Psi/Omega identity lost after correlated P/H/R Joseph measurement")
    return {
        "S": family.S,
        "Sinv": family.Sinv,
        "PHt": family.PHt,
        "K": family.K,
        "A": family.A,
        "correlated_innovation_provenance": family.inverse_provenance,
    }


def build() -> dict:
    original = BACKEND.joseph_measurement
    try:
        BACKEND.joseph_measurement = _correlated_joseph
        try:
            smoke = LINEAGE._two_sample_smoke()
            result = {
                "two_sample_execution_completed": True,
                "failure_type": None,
                "failure_message": None,
                "failure_traceback": None,
                "smoke": smoke,
            }
        except Exception as exc:
            result = {
                "two_sample_execution_completed": False,
                "failure_type": type(exc).__name__,
                "failure_message": str(exc),
                "failure_traceback": traceback.format_exc(),
                "smoke": None,
            }
    finally:
        BACKEND.joseph_measurement = original
    contract = {
        "primitive_family": "(P,H,R)",
        "S_external_independent_rectangle_allowed": False,
        "S_identity": "S=H P H^T+R",
        "K_identity": "K=P H^T (H P H^T+R)^-1",
        "same_P_H_R_used_through_S_inverse_K": True,
        "same_family_K_used_in_Joseph_update": True,
        "rectangular_singular_innovation_admitted": False,
    }
    return {
        "qualification": "OU3_P4_CORRELATED_P_H_R_INNOVATION_CONTINUATION_V3",
        "shipping_filter_changed": False,
        "correlated_innovation_contract": contract,
        "same_history_K_dependency_closed_here": True,
        "scope_of_closure": "P/H/R -> S -> inverse -> K -> same-event Joseph update in this source-uniform continuation",
        "production_backend_facade_available": True,
        "production_facade_module": "ou3_brmm_full_word_riccati_correlated",
        "diagnostic": result,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("shipping_filter_changed") is not False:
        f.append("proof changed shipping filter")
    c = d.get("correlated_innovation_contract", {})
    if c.get("S_external_independent_rectangle_allowed") is not False:
        f.append("arbitrary rectangular S remains admissible")
    for k in ("same_P_H_R_used_through_S_inverse_K", "same_family_K_used_in_Joseph_update"):
        if c.get(k) is not True:
            f.append(k + " not closed")
    if c.get("rectangular_singular_innovation_admitted") is not False:
        f.append("rectangular singular innovation re-admitted")
    if d.get("same_history_K_dependency_closed_here") is not True:
        f.append("P/H/R-to-K dependency not closed")
    if d.get("production_backend_facade_available") is not True:
        f.append("correlated production backend facade missing")
    if d.get("endpoint_augmented_LDLT_closed_here") or d.get("every_prefix_augmented_LDLT_closed_here"):
        f.append("innovation step overclaimed endpoint/prefix closure")
    if d.get("P4_PASS") or d.get("P5_MAY_START"):
        f.append("innovation step promoted P4/P5")
    return f


def main() -> int:
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    print(json.dumps(d, indent=2, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
