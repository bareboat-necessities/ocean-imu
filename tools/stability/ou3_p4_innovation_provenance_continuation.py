#!/usr/bin/env python3
"""Non-promoting continuation past the rectangular innovation-pivot obstruction.

This diagnostic executes the existing two-sample source-uniform H18 smoke path
with only the measurement inverse operation replaced by the already-validated
3x3 PSD-plus-R enclosure.  It does NOT change the canonical Riccati backend and
cannot promote P4/P5.  Its purpose is falsification: determine what obstruction
appears immediately after the false singular member of the entrywise innovation
box is removed.
"""
from __future__ import annotations

import json

from ou3_interval import matrix_add, matrix_identity, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_brmm_full_word_riccati_backend as BACKEND
import ou3_innovation_psd_plus_R_inverse as INNOV
import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE


def _provenance_joseph(state, H, R):
    n = state.dimension
    rows = len(H)
    cols = len(H[0]) if rows else 0
    if cols != n or rows == 0 or len(R) != rows or any(len(row) != rows for row in R):
        raise ValueError("measurement H/R dimension mismatch")
    Ht = matrix_transpose(H)
    PHt = matrix_mul(state.P, Ht)
    S = matrix_symmetric_hull(matrix_add(matrix_mul(H, PHt), R))
    if rows == 3:
        Sinv, inverse_proof = INNOV.innovation_inverse_psd_plus_R_3x3(S, R)
        inverse_mode = "PSD_PLUS_R_3X3_DIAGNOSTIC"
    else:
        Sinv = matrix_inverse_gauss_jordan(S)
        inverse_proof = None
        inverse_mode = "GENERIC_INTERVAL"
    K = matrix_mul(PHt, Sinv)
    A = matrix_sub(matrix_identity(n), matrix_mul(K, H))
    At = matrix_transpose(A)
    KRKt = matrix_mul(matrix_mul(K, R), matrix_transpose(K))
    state.P = matrix_symmetric_hull(matrix_add(matrix_mul(matrix_mul(A, state.P), At), KRKt))
    state.Psi = matrix_mul(A, state.Psi)
    state.Omega = matrix_symmetric_hull(matrix_add(matrix_mul(matrix_mul(A, state.Omega), At), KRKt))
    state.events += 1
    state.measurements += 1
    if not BACKEND.decomposition_identity_enclosed(state):
        raise RuntimeError("P/Psi/Omega identity lost after provenance-aware diagnostic Joseph measurement")
    return {"S": S, "K": K, "A": A, "inverse_mode": inverse_mode, "inverse_proof": inverse_proof}


def build() -> dict:
    original = BACKEND.joseph_measurement
    try:
        BACKEND.joseph_measurement = _provenance_joseph
        try:
            smoke = LINEAGE._two_sample_smoke()
            result = {
                "two_sample_execution_completed": True,
                "failure_type": None,
                "failure_message": None,
                "smoke": smoke,
            }
        except Exception as exc:  # diagnostic records the next fail-closed obstruction
            result = {
                "two_sample_execution_completed": False,
                "failure_type": type(exc).__name__,
                "failure_message": str(exc),
                "smoke": None,
            }
    finally:
        BACKEND.joseph_measurement = original
    return {
        "qualification": "OU3_P4_INNOVATION_PROVENANCE_CONTINUATION_DIAGNOSTIC_V1",
        "canonical_backend_changed": False,
        "shipping_filter_changed": False,
        "rectangular_singular_innovation_admitted": False,
        "same_history_K_dependency_closed_here": False,
        "diagnostic": result,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("canonical_backend_changed") is not False:
        f.append("diagnostic altered canonical backend")
    if d.get("shipping_filter_changed") is not False:
        f.append("diagnostic altered shipping filter")
    if d.get("rectangular_singular_innovation_admitted") is not False:
        f.append("diagnostic re-admitted rectangular singular innovation")
    if d.get("same_history_K_dependency_closed_here") is not False:
        f.append("diagnostic overclaimed same-history K closure")
    if d.get("P4_PASS") or d.get("P5_MAY_START"):
        f.append("diagnostic promoted P4/P5")
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
