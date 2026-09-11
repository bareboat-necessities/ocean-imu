#!/usr/bin/env python3
"""Backward sequential Schur search on the actual H18 post-prediction word.

This is the first certificate search that consumes the corrected sequential
source/finite-precision Schur operator on the exact synchronized shipping word.
It preserves the one physical acceleration endpoint across the word and uses
the coupled a0/a1/J012 sectors at the (single) following prediction.

The search is deliberately fail-closed.  A successful point search is only a
candidate for the source-uniform interval certificate; it does not promote P4.
A failed search reports whether the obstruction is local eliminable-input
positivity or the final entrance-form comparison.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_sub, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_p4_h18_source_indexed_prefix_transport as HTR
import ou3_p4_h18_post_prediction_word as WORD
import ou3_p4_h18_sequential_source_fp_schur as SCHUR
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA = 1
QUALIFICATION = "OU3_P4_H18_SEQUENTIAL_POST_PREDICTION_WORD_SEARCH_V1"
P3_DELTA = 1.0e-18
NP = SCHUR.NP
NX = SCHUR.NX


def I(x: float) -> Interval:
    return Interval.point(float(x))


def zero(n: int):
    return [[I(0.0) for _ in range(n)] for _ in range(n)]


def embed_information(J):
    if len(J) != NX or any(len(r) != NX for r in J):
        raise ValueError("H18 information matrix required")
    Q = zero(NP)
    for i in range(NX):
        for j in range(NX):
            Q[i][j] = J[i][j]
    return Q


def _event_cell_by_token(samples):
    out = {}
    for s in samples:
        for c in s.cells:
            if c.source_token in out:
                raise RuntimeError("duplicate literal event token")
            out[c.source_token] = c
    return out


def _local_delta(kind: str, fp: dict) -> float:
    m = fp["modes"]["H18"]
    dp = float(m["time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction"])
    dm = float(m["explicit_measurement_reset_state_roundoff_norm_upper_per_event"])
    if kind == "prediction":
        return dp
    if kind == "source_coordinate_rebase":
        return 0.0
    return max(dp, dm) if kind == "aw_floor" else dm


def backward_required_form(samples, *, gamma_s: float, gamma_n: float,
                           multipliers: tuple[float, float, float]):
    if gamma_s < 0 or gamma_n <= 0 or any(x < 0 for x in multipliers):
        raise ValueError("nonnegative source and positive finite-precision supplies required")
    word = WORD.build_word(samples)
    fp = FPISS.build()
    ff = FPISS.validate(fp)
    if ff or not fp["additive_ISS_channel_complete_for_conditional_P4"]:
        raise RuntimeError("conditional binary32 ISS prerequisite open")
    by_token = _event_cell_by_token(samples)
    J_terminal = matrix_inverse_gauss_jordan(word["terminal_P"])
    Qnext = embed_information(J_terminal)
    local = []

    for reverse_ordinal, event in enumerate(reversed(word["events"]), start=1):
        kind = event["kind"]
        delta = _local_delta(kind, fp)
        if kind == "prediction":
            cell = by_token.get(event["token"])
            if cell is None:
                raise RuntimeError("terminal prediction token detached from SourceCoverCell")
            h = float(cell.dt_s.hi)
            if not (cell.dt_s.lo == cell.dt_s.hi and h > 0):
                raise RuntimeError("prediction sample time must be one retained point/cell for this search")
            T = SCHUR.prediction_transition(event["C"], cell.tau_applied_s, cell.dt_s, delta)
            sectors_named = SCHUR.prediction_source_sectors(len(T[0]), h=h)
            sectors = [P for _, P in sectors_named]
            prediction = True
        else:
            T = SCHUR.ordinary_transition(event["C"], delta)
            sectors = []
            prediction = False
        try:
            K = SCHUR.local_master(
                Qnext, T, gamma_s=gamma_s, gamma_n=gamma_n,
                sectors=sectors,
                multipliers=multipliers if prediction else (),
                prediction=prediction,
            )
            Qprev, _, piv = SCHUR.schur_required_state_form(K)
        except ValueError as exc:
            return {
                "closed": False,
                "failure_class": "C",
                "failure_stage": "local_Kuu_positive_definiteness",
                "kind": kind,
                "reverse_ordinal": reverse_ordinal,
                "reason": str(exc),
                "local_records": local,
            }
        local.append({
            "kind": kind,
            "reverse_ordinal": reverse_ordinal,
            "Kuu_pivot_lower_min": min(piv),
        })
        Qnext = Qprev

    J_entrance = matrix_inverse_gauss_jordan(word["entrance_P"])
    return {
        "closed": True,
        "required_entrance_form": Qnext,
        "entrance_information": J_entrance,
        "word": word,
        "local_records": local,
    }


def entrance_margin(required, J0, rho: float):
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must lie strictly inside (0,1)")
    desired = zero(NP)
    for i in range(NX):
        for j in range(NX):
            desired[i][j] = I(rho) * J0[i][j]
    margin = matrix_symmetric_hull(matrix_sub(desired, required))
    ok, piv = symmetric_positive_definite_ldlt(margin)
    return ok, [float(x.lo) for x in piv]


def candidate_grid():
    # Deterministic, deliberately broad.  These are search values only; a
    # successful tuple still needs source-uniform outward subdivision.
    rhos = (0.90, 0.95, 0.98, 0.99, 0.995, 0.999, 0.9999)
    gammas_n = tuple(10.0 ** p for p in range(-6, 25, 3))
    gammas_s = tuple(10.0 ** p for p in range(-6, 25, 3))
    lambdas = tuple(10.0 ** p for p in range(-6, 25, 3))
    return rhos, gammas_s, gammas_n, lambdas


def search(samples):
    rhos, gs_values, gn_values, lam_values = candidate_grid()
    attempts = 0
    first_local_failure = None
    best = None
    for lam in lam_values:
        multipliers = (lam, lam, lam)
        for gn in gn_values:
            # gamma_s only affects the prediction homogeneous source scale; try
            # smallest first so the reported candidate remains informative.
            for gs in gs_values:
                rec = backward_required_form(samples, gamma_s=gs, gamma_n=gn,
                                             multipliers=multipliers)
                attempts += 1
                if not rec["closed"]:
                    if first_local_failure is None:
                        first_local_failure = {
                            k: v for k, v in rec.items()
                            if k not in ("local_records", "required_entrance_form", "word")
                        }
                    continue
                for rho in rhos:
                    ok, piv = entrance_margin(rec["required_entrance_form"], rec["entrance_information"], rho)
                    margin = min(piv) if piv else -math.inf
                    if best is None or margin > best["pivot_lower_min"]:
                        best = {
                            "rho": rho, "gamma_s": gs, "gamma_n": gn,
                            "multipliers": list(multipliers),
                            "pivot_lower_min": margin,
                            "entrance_LDLT_closed": ok,
                        }
                    if ok:
                        return {
                            "candidate_found": True,
                            "attempts": attempts,
                            "candidate": best,
                            "first_local_failure": first_local_failure,
                            "event_count": len(rec["word"]["events"]),
                            "local_Kuu_pivot_lower_min": min(x["Kuu_pivot_lower_min"] for x in rec["local_records"]),
                        }
    return {
        "candidate_found": False,
        "attempts": attempts,
        "candidate": best,
        "first_local_failure": first_local_failure,
    }


def build():
    schur = SCHUR.build(); sf = SCHUR.validate(schur)
    word = WORD.build(); wf = WORD.validate(word)
    fp = FPISS.build(); ff = FPISS.validate(fp)
    if sf or wf or ff:
        raise RuntimeError(f"search prerequisites failed Schur={sf} word={wf} fp={ff}")
    result = search(HTR.smoke_objects())
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "P3_delta": P3_DELTA,
        "actual_post_prediction_shipping_word_consumed": True,
        "actual_entrance_and_terminal_P_inverse_consumed": True,
        "same_history_physical_acceleration_endpoint_carried": True,
        "coupled_a0_a1_J012_sectors_consumed_at_prediction": True,
        "conditional_binary32_event_channel_consumed": True,
        "nonnegative_S_procedure_multiplier_search_executed": True,
        "strict_endpoint_rho_search_below_one_executed": True,
        "search": result,
        "point_candidate_is_not_source_uniform_certificate": True,
        "source_uniform_outward_subdivision_closed_here": False,
        "production_endpoint_augmented_LDLT_closed_here": False,
        "production_every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "if the sequential point search finds a candidate, freeze its topology only and rerun with source-dependent interval cells/subdivision; "
            "if it fails, use the recorded local Kuu or entrance margin to refine multipliers/storage without interpreting dependency loss as instability"
        ),
    }


def validate(d):
    f = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    for k in (
        "actual_post_prediction_shipping_word_consumed",
        "actual_entrance_and_terminal_P_inverse_consumed",
        "same_history_physical_acceleration_endpoint_carried",
        "coupled_a0_a1_J012_sectors_consumed_at_prediction",
        "conditional_binary32_event_channel_consumed",
        "nonnegative_S_procedure_multiplier_search_executed",
        "strict_endpoint_rho_search_below_one_executed",
        "point_candidate_is_not_source_uniform_certificate",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "source_uniform_outward_subdivision_closed_here",
        "production_endpoint_augmented_LDLT_closed_here",
        "production_every_prefix_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    s = d.get("search", {})
    if not isinstance(s.get("attempts"), int) or s["attempts"] <= 0:
        f.append("sequential search did not execute")
    return list(dict.fromkeys(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    d = build(); f = validate(d)
    d["validation_pass"] = not f; d["validation_failures"] = f
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(d, indent=2, sort_keys=True, default=str) + "\n")
    s = d["search"]
    print(json.dumps({
        "candidate": s.get("candidate_found"),
        "attempts": s.get("attempts"),
        "best": s.get("candidate"),
        "P4": d["P4_PASS"], "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
