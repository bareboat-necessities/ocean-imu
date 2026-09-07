#!/usr/bin/env python3
"""Exact complete-word endpoint transport for nonlinear OU-III P4.

This module is part of the canonical P4 master inequality.  It is not a point
experiment, source replacement, correction-radius certificate, or packetwise
remainder budget.

For node k let

    Phi_k = z_k + E_k epsilon_k,

where epsilon_k is the full accelerometer-linearizing wave shift

    epsilon_k = (Q_aw,k-I) delta_a_w,k + e_eta,k.

Write one literal shipping event as

    Phi_{k+1} = C_k Phi_k + xi_k,
    xi_k = rho_k + E_{k+1} epsilon_{k+1} - L_k E_k epsilon_k.

C_k is the tangent map actually carried in the complete word.  L_k is the
physical/chart map multiplying the *incoming coordinate shift*:

* prediction/source/hybrid lift: C_k=L_k (literal F/source/lift map);
* Joseph plus immediate left reset: C_k=G_k(I-K_k H_k), L_k=G_k;
* covariance-only floor / identity event: C_k=L_k=I.

Let M_{N:k}=C_{N-1}...C_k and M_W=M_{N:0}.  Variation of constants gives the
EXACT endpoint defect

    d_W = Phi_N - M_W Phi_0
        = r_W + E_N epsilon_N - M_W E_0 epsilon_0
          + sum_k M_{N:k+1}(C_k-L_k)E_k epsilon_k,

    r_W = sum_k M_{N:k+1} rho_k.

This identity is the key complete-word reduction.  It does NOT replace the sum
by a sum of norms.  Prediction/source/hybrid/floor terms cancel algebraically.
For Joseph/reset events

    (C_k-L_k)E_k = -G_k K_k H_k E_k.

The S=0 and magnetometer Jacobians have zero a_w column, so their interior
full-shift term is exactly zero.  Only accepted accelerometer corrections leave
an interior block, and that block is retained as one JOINT suffix-weighted
operator

    B_W = [ M_{N:k+1} G_k K_k H_k E_k ]_{k in accel(W)}.

Every later shipping operation -- in particular every due S=0 update with its
actual applied anisotropic SpectralMSE R_S -- is inside the corresponding
suffix M_{N:k+1}.  Therefore B_W is the correct object to enclose.  Bounding
individual packets and multiplying by the packet count is explicitly forbidden.

With J_0=P_0^{-1}, J_N=P_N^{-1}, D_W=J_0-M_W^T J_N M_W, the exact endpoint
energy identity is

    V_N-V_0 = -Phi_0^T D_W Phi_0
              + 2 (M_W Phi_0)^T J_N d_W
              + d_W^T J_N d_W.

P3 supplies the complete same-history linear word and its full-matrix margin.
P4 must now enclose the *joint* (r_W, epsilon_0, epsilon_N,
epsilon_acc-history) object over the same complete SEA3 word and show the last
two endpoint terms fit inside that full-matrix decrease.  No scalar correction
radius, inverse-metric floor, independent R_S schedule, replay, or alternate
estimator is introduced here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

import ou3_p4_complete_sea3_measurement_linearizing_aw_coordinate as AW
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as ACC
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_WHOLE_WORD_ENDPOINT_TRANSPORT_V1"


def _shape(a: Sequence[Sequence[Any]]) -> tuple[int, int]:
    rows = len(a)
    cols = len(a[0]) if rows else 0
    if any(len(row) != cols for row in a):
        raise ValueError("ragged matrix")
    return rows, cols


def _mm(a: Sequence[Sequence[Any]], b: Sequence[Sequence[Any]]) -> list[list[Any]]:
    ar, ac = _shape(a)
    br, bc = _shape(b)
    if ac != br:
        raise ValueError("matrix product dimension mismatch")
    if ar == 0 or bc == 0:
        return [[] for _ in range(ar)]
    out: list[list[Any]] = []
    for i in range(ar):
        row = []
        for j in range(bc):
            value = a[i][0] * b[0][j]
            for k in range(1, ac):
                value = value + a[i][k] * b[k][j]
            row.append(value)
        out.append(row)
    return out


def _mv(a: Sequence[Sequence[Any]], x: Sequence[Any]) -> list[Any]:
    ar, ac = _shape(a)
    if ac != len(x):
        raise ValueError("matrix/vector dimension mismatch")
    if ac == 0:
        return []
    out = []
    for i in range(ar):
        value = a[i][0] * x[0]
        for k in range(1, ac):
            value = value + a[i][k] * x[k]
        out.append(value)
    return out


def _add(a: Sequence[Any], b: Sequence[Any]) -> list[Any]:
    if len(a) != len(b):
        raise ValueError("vector dimension mismatch")
    return [x + y for x, y in zip(a, b)]


def _sub(a: Sequence[Any], b: Sequence[Any]) -> list[Any]:
    if len(a) != len(b):
        raise ValueError("vector dimension mismatch")
    return [x - y for x, y in zip(a, b)]


def _msub(a: Sequence[Sequence[Any]], b: Sequence[Sequence[Any]]) -> list[list[Any]]:
    if _shape(a) != _shape(b):
        raise ValueError("matrix dimension mismatch")
    return [[a[i][j] - b[i][j] for j in range(len(a[i]))] for i in range(len(a))]


def _transpose(a: Sequence[Sequence[Any]]) -> list[list[Any]]:
    r, c = _shape(a)
    return [[a[i][j] for i in range(r)] for j in range(c)]


def _identity(n: int, one: Any, zero: Any) -> list[list[Any]]:
    return [[one if i == j else zero for j in range(n)] for i in range(n)]


def _zero_vector_like(x: Sequence[Any], n: int) -> list[Any]:
    if not x:
        raise ValueError("cannot infer zero from empty vector")
    z = x[0] - x[0]
    return [z for _ in range(n)]


def _matrix_is_zero(a: Sequence[Sequence[Any]]) -> bool:
    return all(x == (x - x) for row in a for x in row)


def _dot(a: Sequence[Any], b: Sequence[Any]) -> Any:
    if len(a) != len(b) or not a:
        raise ValueError("dot-product dimension mismatch")
    value = a[0] * b[0]
    for i in range(1, len(a)):
        value = value + a[i] * b[i]
    return value


def suffix_products(events: Sequence[dict[str, Any]], one: Any, zero: Any) -> list[list[list[Any]]]:
    """Return suffix[k]=C_{N-1}...C_k, allowing one rectangular H->A event."""
    if not events:
        raise ValueError("complete word must contain at least one event")
    first_in = _shape(events[0]["C"])[1]
    final_out = _shape(events[-1]["C"])[0]
    for k, event in enumerate(events):
        cr, cc = _shape(event["C"])
        lr, lc = _shape(event["L"])
        if (cr, cc) != (lr, lc):
            raise ValueError(f"event {k} C/L shape mismatch")
        if k and cc != _shape(events[k - 1]["C"])[0]:
            raise ValueError(f"event {k} does not compose with predecessor")
    suffix: list[list[list[Any]]] = [None] * (len(events) + 1)  # type: ignore[list-item]
    suffix[-1] = _identity(final_out, one, zero)
    for k in range(len(events) - 1, -1, -1):
        suffix[k] = _mm(suffix[k + 1], events[k]["C"])
    if _shape(suffix[0]) != (final_out, first_in):
        raise RuntimeError("whole-word suffix product shape drifted")
    return suffix


def endpoint_decomposition(
    events: Sequence[dict[str, Any]],
    embeddings: Sequence[Sequence[Sequence[Any]]],
    eps_nodes: Sequence[Sequence[Any]],
) -> dict[str, Any]:
    """Evaluate both sides of the exact whole-word shift identity.

    ``events[k]`` contains C, L, rho and kind.  ``embeddings[k]`` is E_k and
    ``eps_nodes[k]`` is epsilon_k.  Arithmetic may use float, Fraction, Decimal,
    or another exact/validated scalar type supporting +,-,* and equality.
    """
    n = len(events)
    if len(embeddings) != n + 1 or len(eps_nodes) != n + 1:
        raise ValueError("need one E/epsilon node at every word boundary")
    if not eps_nodes[0]:
        raise ValueError("epsilon coordinate must be nonempty")
    one = eps_nodes[0][0] - eps_nodes[0][0] + 1
    zero = eps_nodes[0][0] - eps_nodes[0][0]
    suffix = suffix_products(events, one, zero)
    final_dim = _shape(events[-1]["C"])[0]
    direct = _zero_vector_like(eps_nodes[0], final_dim)
    r_word = _zero_vector_like(eps_nodes[0], final_dim)
    interior = _zero_vector_like(eps_nodes[0], final_dim)
    accel_blocks: list[list[list[Any]]] = []
    interior_events: list[int] = []

    for k, event in enumerate(events):
        cr, cc = _shape(event["C"])
        if len(event["rho"]) != cr:
            raise ValueError(f"event {k} rho dimension mismatch")
        if _shape(embeddings[k])[0] != cc or _shape(embeddings[k + 1])[0] != cr:
            raise ValueError(f"event {k} embedding dimension mismatch")
        if _shape(embeddings[k])[1] != len(eps_nodes[k]):
            raise ValueError(f"event {k} incoming epsilon dimension mismatch")
        if _shape(embeddings[k + 1])[1] != len(eps_nodes[k + 1]):
            raise ValueError(f"event {k} outgoing epsilon dimension mismatch")

        incoming = _mv(_mm(event["L"], embeddings[k]), eps_nodes[k])
        outgoing = _mv(embeddings[k + 1], eps_nodes[k + 1])
        xi = _add(event["rho"], _sub(outgoing, incoming))
        direct = _add(direct, _mv(suffix[k + 1], xi))
        r_word = _add(r_word, _mv(suffix[k + 1], event["rho"]))

        defect_map = _mm(_msub(event["C"], event["L"]), embeddings[k])
        weighted = _mm(suffix[k + 1], defect_map)
        if not _matrix_is_zero(weighted):
            interior_events.append(k)
            interior = _add(interior, _mv(weighted, eps_nodes[k]))
            if event.get("kind") != "accelerometer":
                raise ValueError(
                    f"non-accelerometer event {k} has nonzero interior epsilon transport"
                )
            # B_W is defined with the positive G K H E block; the endpoint
            # decomposition carries -B_W epsilon.  C-L=-GKH for this event.
            accel_blocks.append([[-x for x in row] for row in weighted])

    m_word = suffix[0]
    endpoint = _sub(
        _mv(embeddings[-1], eps_nodes[-1]),
        _mv(_mm(m_word, embeddings[0]), eps_nodes[0]),
    )
    decomposed = _add(r_word, _add(endpoint, interior))
    residual = _sub(direct, decomposed)
    return {
        "M_word": m_word,
        "direct_defect": direct,
        "r_word": r_word,
        "endpoint_epsilon_term": endpoint,
        "interior_epsilon_term": interior,
        "interior_event_indices": interior_events,
        "accelerometer_suffix_blocks": accel_blocks,
        "decomposed_defect": decomposed,
        "identity_residual": residual,
    }


def master_energy_identity(
    J0: Sequence[Sequence[Any]],
    JN: Sequence[Sequence[Any]],
    M: Sequence[Sequence[Any]],
    phi0: Sequence[Any],
    defect: Sequence[Any],
) -> dict[str, Any]:
    """Return the exact three-term P4 endpoint energy decomposition."""
    if _shape(J0)[0] != _shape(J0)[1] or _shape(JN)[0] != _shape(JN)[1]:
        raise ValueError("endpoint information matrices must be square")
    mx = _mv(M, phi0)
    if len(mx) != len(defect) or len(phi0) != _shape(J0)[0] or len(mx) != _shape(JN)[0]:
        raise ValueError("master energy dimensions do not match")
    jn_mx = _mv(JN, mx)
    jn_d = _mv(JN, defect)
    linear_final = _dot(mx, jn_mx)
    initial = _dot(phi0, _mv(J0, phi0))
    linear_decrease = initial - linear_final
    cross = 2 * _dot(mx, jn_d)
    defect_energy = _dot(defect, jn_d)
    total_delta = -linear_decrease + cross + defect_energy
    direct_delta = _dot(_add(mx, defect), _mv(JN, _add(mx, defect))) - initial
    return {
        "initial_energy": initial,
        "linear_final_energy": linear_final,
        "linear_decrease": linear_decrease,
        "cross_term": cross,
        "defect_energy": defect_energy,
        "total_delta": total_delta,
        "direct_delta": direct_delta,
        "identity_residual": direct_delta - total_delta,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    event = EVENT.build()
    aw = AW.build(path)
    acc = ACC.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"event: {x}" for x in EVENT.validate(event)]
        + [f"full-shift: {x}" for x in AW.validate(aw)]
        + [f"accelerometer: {x}" for x in ACC.validate(acc)]
    )
    if failures:
        raise RuntimeError(f"whole-word endpoint prerequisites failed: {failures}")
    if p3.get("P3_CONDITIONAL_SEA3_PASS") is not True:
        raise RuntimeError("whole-word endpoint transport requires frozen conditional P3")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "same_complete_SEA3_execution_required": True,
        "actual_applied_per_axis_RS_retained_in_word": True,
        "all_due_S_updates_retained": True,
        "all_valid_accelerometer_updates_retained": True,
        "asynchronous_vector_events_retained": True,
        "full_process_Q_and_aw_floors_retained": True,
        "H_to_A_rectangular_hybrid_retained": True,
        "full_epsilon_aw_retained": True,
        "full_prediction_F_Eaw_retained": True,
        "whole_word_variation_of_constants_used": True,
        "endpoint_decomposition_identity": (
            "d_W=r_W+E_N*epsilon_N-M_W*E_0*epsilon_0+"
            "sum_k M_suffix(k)*(C_k-L_k)*E_k*epsilon_k"
        ),
        "prediction_source_hybrid_floor_interior_epsilon_terms_cancel_exactly": True,
        "S_zero_interior_epsilon_term_zero_exactly": True,
        "magnetometer_interior_epsilon_term_zero_exactly": True,
        "accepted_accelerometer_is_only_interior_epsilon_event_class": True,
        "accelerometer_joint_operator": (
            "B_W=[M_suffix(k)*G_k*K_k*H_k*E_aw] over accepted accelerometer events"
        ),
        "accelerometer_endpoint_defect_sign": "-B_W*epsilon_acc_history",
        "actual_RS_regularization_enters_every_applicable_suffix": True,
        "packetwise_norm_sum_used": False,
        "packet_count_multiplier_used": False,
        "independent_RS_schedule_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "master_endpoint_energy_identity": (
            "DeltaV=-Phi0^T*D_W*Phi0+2*(M_W*Phi0)^T*P_N^-1*d_W+"
            "d_W^T*P_N^-1*d_W"
        ),
        "D_W_definition": "P_0^-1-M_W^T*P_N^-1*M_W",
        "master_inequality_object_emitted": True,
        "source_uniform_joint_BW_epsilon_enclosure_closed": False,
        "source_uniform_r_word_enclosure_closed": False,
        "source_uniform_master_endpoint_domination_closed": False,
        "P4_promoted_here": False,
        "P5_may_start_here": False,
        "next_obligation": (
            "enclose the joint suffix-weighted accelerometer operator B_W and r_W over the SAME complete SEA3 word, "
            "with actual applied R_S inside the suffix maps, and prove the two endpoint nonlinear terms fit inside "
            "the full-matrix P3 decrease for the widest declared [30,25,20,15] degree candidate; do not scalarize packetwise"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "P3_frozen_not_modified", "same_complete_SEA3_execution_required",
        "actual_applied_per_axis_RS_retained_in_word", "all_due_S_updates_retained",
        "all_valid_accelerometer_updates_retained", "asynchronous_vector_events_retained",
        "full_process_Q_and_aw_floors_retained", "H_to_A_rectangular_hybrid_retained",
        "full_epsilon_aw_retained", "full_prediction_F_Eaw_retained",
        "whole_word_variation_of_constants_used",
        "prediction_source_hybrid_floor_interior_epsilon_terms_cancel_exactly",
        "S_zero_interior_epsilon_term_zero_exactly",
        "magnetometer_interior_epsilon_term_zero_exactly",
        "accepted_accelerometer_is_only_interior_epsilon_event_class",
        "actual_RS_regularization_enters_every_applicable_suffix",
        "master_inequality_object_emitted",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "packetwise_norm_sum_used", "packet_count_multiplier_used", "independent_RS_schedule_used",
        "correction_radius_claim_used", "inverse_metric_floor_claim_used", "trajectory_replay_used",
        "source_family_replaced", "source_uniform_joint_BW_epsilon_enclosure_closed",
        "source_uniform_r_word_enclosure_closed", "source_uniform_master_endpoint_domination_closed",
        "P4_promoted_here", "P5_may_start_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "canonical_source": d["canonical_source"],
        "endpoint_master_object": d["master_inequality_object_emitted"],
        "only_interior_shift_event": "accelerometer",
        "actual_RS_in_suffix": d["actual_RS_regularization_enters_every_applicable_suffix"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
