#!/usr/bin/env python3
"""Correlation-preserving facade for the canonical OU-III Riccati backend.

All non-measurement operations delegate to ou3_brmm_full_word_riccati_backend.
Measurement operations are replaced by one primitive that carries the same
(P,H,R) family through S=HPH^T+R, S^-1 and K.  S is never accepted as an
independent caller-supplied rectangle.

This facade is used by the literal H18/A21 word assembler.  It is intentionally
small so the algebraic P/Psi/Omega backend remains unchanged while the source
proof stops breaking the innovation dependency at the Joseph boundary.
"""
from __future__ import annotations

import ou3_brmm_full_word_riccati_backend as OLD
import ou3_correlated_innovation_family as CORR
from ou3_interval import matrix_add, matrix_mul, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull

USEFUL_GATE = OLD.USEFUL_GATE
JointWordState = OLD.JointWordState
PriorFreeBatchState = OLD.PriorFreeBatchState

# Explicitly re-export the ordinary backend operations used by the word assembler.
initialize = OLD.initialize
initialize_prior_free = OLD.initialize_prior_free
predict = OLD.predict
prior_free_predict = OLD.prior_free_predict
add_psd_floor = OLD.add_psd_floor
prior_free_add_prior_independent_psd_floor = OLD.prior_free_add_prior_independent_psd_floor
reconstruct_joint_from_prior_free = OLD.reconstruct_joint_from_prior_free
uniform_prior_completion_matrix = OLD.uniform_prior_completion_matrix
certify_uniform_prior_contraction = OLD.certify_uniform_prior_contraction
certify_contraction = OLD.certify_contraction
decomposition_identity_enclosed = OLD.decomposition_identity_enclosed
contraction_matrix = OLD.contraction_matrix
prediction_contraction_image = OLD.prediction_contraction_image
joseph_contraction_image = OLD.joseph_contraction_image
floor_contraction_image = OLD.floor_contraction_image
contraction_identity_enclosed = OLD.contraction_identity_enclosed
shipping_source_parity = OLD.shipping_source_parity
contraction_preservation_identities = OLD.contraction_preservation_identities
prior_free_batch_identities = OLD.prior_free_batch_identities


def joseph_measurement(state: JointWordState, H, R) -> dict:
    """Joseph update from one correlated (P,H,R) innovation family."""
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
    if not OLD.decomposition_identity_enclosed(state):
        raise RuntimeError("P/Psi/Omega identity lost after correlated Joseph measurement")
    return {
        "S": family.S,
        "Sinv": family.Sinv,
        "K": family.K,
        "A": family.A,
        "PHt": family.PHt,
        "correlated_innovation_provenance": family.inverse_provenance,
    }


def prior_free_measurement(state: PriorFreeBatchState, H, R) -> dict:
    """Prior-free Joseph information update with correlated (Qc,H,R) family."""
    family = CORR.build(state.Qc, H, R)
    J = matrix_mul(H, state.T)
    state.D = matrix_symmetric_hull(
        matrix_add(
            state.D,
            matrix_mul(matrix_mul(matrix_transpose(J), family.Sinv), J),
        )
    )
    K0, A0 = family.K, family.A
    A0t = matrix_transpose(A0)
    K0RK0t = matrix_mul(matrix_mul(K0, family.R), matrix_transpose(K0))
    state.T = matrix_mul(A0, state.T)
    state.Qc = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A0, state.Qc), A0t), K0RK0t)
    )
    state.events += 1
    state.measurements += 1
    return {
        "S0": family.S,
        "S0inv": family.Sinv,
        "K0": family.K,
        "A0": family.A,
        "J": J,
        "QcHt": family.PHt,
        "correlated_innovation_provenance": family.inverse_provenance,
    }


def correlated_innovation_contract() -> dict:
    return {
        "primitive_family": "(P,H,R)",
        "S_external_independent_rectangle_allowed": False,
        "same_P_H_R_used_through_S_inverse_K": True,
        "formula": "K=P H^T (H P H^T+R)^-1",
        "three_component_invertibility_uses_PSD_plus_R": True,
        "shipping_filter_changed": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate_backend() -> list[str]:
    f = list(OLD.validate_backend())
    reg = CORR.validate_point_regression()
    if reg.get("same_P_H_R_used_for_PHt_S_Sinv_K") is not True:
        f.append("correlated innovation family lost P/H/R ancestry")
    if reg.get("S_is_external_independent_argument") is not False:
        f.append("innovation S can still be injected independently")
    if reg.get("rectangular_S_family_can_be_selected_independently") is not False:
        f.append("rectangular S family can still be selected independently")
    return list(dict.fromkeys(f))


if __name__ == "__main__":
    import json
    f = validate_backend()
    print(json.dumps({"contract": correlated_innovation_contract(), "failures": f}, indent=2, sort_keys=True))
    raise SystemExit(0 if not f else 2)
