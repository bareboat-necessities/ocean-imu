#!/usr/bin/env python3
"""Audit retained full-word point witnesses against conditional BIAS0/1/2.

No new direction, scale, source, metric or estimator is searched. The inputs
are the archived directions/scales and content-addressed physical payloads.
The existing finite-value evaluator and reset-normalized covariance are reused.
Source membership and finite reset/storage attachment are separate obligations:
a zero true-bias root satisfies GM algebra but does not bind a replay payload to
the complete BRMM joint source. Hardware qualification is not a prerequisite
for this conditional experiment. Numerical point ratios are not certificates.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import tarfile
from pathlib import Path

import numpy as np

import ou3_p4_physical_finite_map_feasibility as BASE
import ou3_p4_physical_finite_map_feasibility_fast as FAST
import ou3_mems_bias_contract as BIAS
import ou3_p4_retained_word_attachment as ATTACH

HERE = Path(__file__).resolve().parent
DEFAULT_WITNESSES = HERE / "fixtures/ou3_p4_retained_bias_witnesses.json"


def unpack_retained_payloads(archive: Path, prefix: Path, witness: dict):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != witness["payload_archive_sha256"]:
        raise ValueError("retained payload archive identity changed")
    expected = {
        f"ou3_p4_physical_payload{suffix}.bin": Path(f"{prefix}{suffix}.bin")
        for mode in ("H18", "A21") for suffix in (f"_{mode}", f"_canonical_{mode}")
    }
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:xz") as stream:
        members = stream.getmembers()
        if len(members) != 4 or {m.name for m in members} != set(expected) or not all(m.isfile() for m in members):
            raise ValueError("retained archive must contain exactly the four regular payload files")
        for member in members:
            with stream.extractfile(member) as source:
                expected[member.name].write_bytes(source.read())


TRACE_FIELDS = (
    "mode", "scale", "event_index", "event", "time", "elapsed_prediction_s",
    "estimated_bias_norm", "bias_error_norm", "error_bx", "error_by", "error_bz",
    "g_x", "g_y", "g_z", "y_x", "y_y", "y_z", "projection_defect_norm",
    "energy_ratio",
)


def separation_ratios(a_energy: float, d_energy: float, cross: float,
                      y_energy: float, initial_energy: float) -> dict:
    """BIAS2 point evidence with Xi=V0, retaining the dense history cross term.

These quotients are evaluated on one nonlinear trajectory. Taking their minimum
over the archive does not establish a bound on the source/error graph.
"""
    values = (a_energy, d_energy, cross, y_energy, initial_energy)
    if not all(math.isfinite(v) for v in values) or initial_energy <= 0:
        raise ValueError("finite energies and positive initial storage required")
    total = a_energy + d_energy
    if min(a_energy, d_energy, y_energy) < 0 or total <= 0:
        raise ValueError("nonnegative channel energies with nonzero total required")
    kappa = 2.0 * abs(cross) / total
    return {
        "Xi": "V0=x0^T J0 x0 in the retained full-state storage",
        "g_weighted_energy": a_energy,
        "corrected_bias_error_weighted_energy": d_energy,
        "weighted_cross_inner_product": cross,
        "residual_weighted_energy": y_energy,
        "energy_identity_residual": y_energy - (total + 2.0 * cross),
        "kappa_point": kappa,
        "alpha_point": total / initial_energy,
        "mu_sufficient_point": (total - 2.0 * abs(cross)) / initial_energy,
        "mu_actual_point": y_energy / initial_energy,
        "positive_sufficient_separation_at_this_point": kappa < 1.0,
        "source_uniform_mu_certified": False,
    }


def case_verdict(bias_interior: bool, legacy_domain: bool,
                 bias_recurrence: bool) -> str:
    if not bias_recurrence:
        return "BIAS1_POINT_RECURRENCE_FAILED"
    if not bias_interior or not legacy_domain:
        return "OUTSIDE_CHECKED_CONDITIONAL_POINT_DOMAIN"
    return "BIAS_COMPATIBLE_FULL_BRMM_ADMISSIBILITY_UNRESOLVED"


def audit_case(payload: dict, linear: dict, kernels: list[dict], domain: dict,
               scale: float, source_input: Path, trace=None) -> dict:
    n = payload["mode_dim"]
    mode = "H18" if n == 18 else "A21"
    state = float(scale) * np.asarray(linear["direction"], dtype=float)
    V0 = float(state @ linear["Q0"] @ state)
    if not V0 > 0.0 or not math.isfinite(V0):
        raise ValueError("retained witness has invalid initial storage")
    radius = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    interior = float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
    hs = BASE.parse_hs(source_input)
    initial_error = state[18:21].copy() if n == 21 else np.zeros(3)
    free_error = initial_error.copy()
    elapsed = 0.0
    max_estimate = 0.0
    max_free_error_defect = 0.0
    max_recurrence_defect = 0.0
    max_correction = 0.0
    max_projection = 0.0
    max_ratio = 1.0
    first_interior_failure = None
    first_domain_failure = None
    projection_counts = {"active": 0, "inactive": 0, "clarke_hull": 0, "not_applicable": 0}
    aa = dd = ad = yy = 0.0
    finite_reset = {"event_count": 0, "max_defect_norm": 0.0,
                    "max_relative_signed_identity_residual": 0.0,
                    "max_absolute_attached_minus_unreset_storage": 0.0,
                    "whole_word_storage_attached": False}

    def retain_state(index, event, time, g, y, proj):
        nonlocal max_estimate, first_interior_failure, first_domain_failure, max_ratio
        eb = state[18:21] if n == 21 else np.zeros(3)
        # The archived homogeneous point evaluator fixes b_true=beta_root=0.
        # Therefore b_hat=-e_b exactly; a new per-event true bias is never chosen.
        bnorm = float(np.linalg.norm(eb))
        max_estimate = max(max_estimate, bnorm)
        if n == 21 and bnorm > interior and first_interior_failure is None:
            first_interior_failure = {
                "event_index": index, "event": event, "time": time,
                "estimated_bias_norm": bnorm, "bound": interior,
            }
        ok, reason = BASE.domain_status(state, n, domain, hs)
        if not ok and first_domain_failure is None:
            first_domain_failure = {"event_index": index, "event": event, "time": time, "reason": reason}
        J = linear["Q0"] if index < 0 else kernels[index]["Qafter"]
        ratio = float(state @ J @ state) / V0
        max_ratio = max(max_ratio, ratio)
        if trace is not None:
            trace.writerow(dict(zip(TRACE_FIELDS, (
                mode, scale, index, event, time, elapsed, bnorm, float(np.linalg.norm(eb)),
                *eb, *g, *y, proj, ratio,
            ))))

    retain_state(-1, "initial", payload["t0"], np.zeros(3), np.zeros(3), 0.0)
    for idx, (ev, kernel) in enumerate(zip(payload["events"], kernels, strict=True)):
        before = state.copy()
        eb = before[18:21] if n == 21 else np.zeros(3)
        g = y = np.zeros(3)
        proj_norm = 0.0
        if ev["type"] == BASE.EV_PRED:
            state = FAST._prediction_value(ev, before, n)
            elapsed += ev["h"]
            if n == 21:
                phi = float(ev["linear_shipping"][18, 18])
                expected = phi * eb
                free_error = phi * free_error
                max_recurrence_defect = max(max_recurrence_defect, float(np.linalg.norm(state[18:21] - expected)))
        elif ev["type"] == BASE.EV_FLOOR:
            state = before.copy()
        else:
            y = FAST._physical_residual(ev, before, n)
            K = kernel["K"]
            if K is None:
                raise ValueError("measurement lost its retained Joseph gain")
            if ev["type"] == BASE.EV_ACC:
                g = y - eb
                # R is the same event covariance used in the retained gain.
                aa += float(g @ np.linalg.solve(ev["R"], g))
                dd += float(eb @ np.linalg.solve(ev["R"], eb))
                ad += float(g @ np.linalg.solve(ev["R"], eb))
                yy += float(y @ np.linalg.solve(ev["R"], y))
            state, branch = FAST._joseph_value(ev, before, n, K, radius)
            # Attach each finite physical reset to the covariance congruence
            # generated by ITS correction. This does not justify reusing the
            # following reset-deleted P/K cells as a shipping word.
            attachment = ATTACH.finite_reset_energy(
                before, state, K, y, linear["path"][idx]["Pafter"])
            finite_reset["event_count"] += 1
            finite_reset["max_defect_norm"] = max(finite_reset["max_defect_norm"], attachment["reset_defect_norm"])
            finite_reset["max_relative_signed_identity_residual"] = max(
                finite_reset["max_relative_signed_identity_residual"], attachment["relative_signed_identity_residual"])
            finite_reset["max_absolute_attached_minus_unreset_storage"] = max(
                finite_reset["max_absolute_attached_minus_unreset_storage"], abs(attachment["attached_minus_unreset_storage"]))
            projection_counts[branch] += 1
            if n == 21:
                correction = K[18:21] @ y
                estimate_corr = -eb + correction
                norm_corr = float(np.linalg.norm(estimate_corr))
                estimate_proj = estimate_corr if norm_corr <= radius else estimate_corr * (radius / norm_corr)
                projection_defect = estimate_corr - estimate_proj
                expected = eb - correction + projection_defect
                max_recurrence_defect = max(max_recurrence_defect, float(np.linalg.norm(state[18:21] - expected)))
                max_correction = max(max_correction, float(np.linalg.norm(correction)))
                proj_norm = float(np.linalg.norm(projection_defect))
                max_projection = max(max_projection, proj_norm)
        if n == 21:
            max_free_error_defect = max(max_free_error_defect, float(np.linalg.norm(state[18:21] - free_error)))
        retain_state(idx, ev["name"], ev["time"], g, y, proj_norm)

    rho = float(state @ linear["QN"] @ state) / V0
    separation = separation_ratios(aa, dd, ad, yy, V0)
    recurrence_ok = max_recurrence_defect <= 1e-12
    interior_ok = first_interior_failure is None
    domain_ok = first_domain_failure is None
    return {
        "scale": scale, "initial_state": (scale * linear["direction"]).tolist(),
        "initial_energy": V0, "rho_endpoint": rho, "distance_to_one": 1.0 - rho,
        "max_prefix_ratio": max_ratio, "endpoint_state": state.tolist(),
        "true_bias_root": [0.0, 0.0, 0.0],
        "true_bias_GM_history_exact_for_any_positive_tau": True,
        "external_bias_forcing_zero": True,
        "initial_and_all_event_bias_interiors_checked": n == 21,
        "max_estimated_bias_norm": max_estimate,
        "estimated_bias_interior_retained": interior_ok,
        "first_estimated_bias_interior_failure": first_interior_failure,
        "legacy_state_domain_retained_including_initial": domain_ok,
        "first_legacy_state_domain_failure": first_domain_failure,
        "max_bias_recurrence_point_residual": max_recurrence_defect,
        "max_measurement_bias_correction": max_correction,
        "max_projection_defect": max_projection,
        "max_corrected_error_minus_free_GM_path": max_free_error_defect,
        "projection_counts": projection_counts,
        "finite_reset_energy_attachment": finite_reset,
        "bias_recurrence_point_check_pass": recurrence_ok,
        "conditional_point_checks_pass": recurrence_ok and interior_ok and domain_ok,
        "BIAS2": separation,
        "positive_point_separation_and_expanding_storage": rho > 1.0 and separation["positive_sufficient_separation_at_this_point"],
        "verdict": case_verdict(interior_ok, domain_ok, recurrence_ok),
        "complete_BRMM_source_and_bias_membership_verified": False,
        "finite_reset_storage_attachment_verified": False,
        "canonical_P4_falsified_here": False,
    }


def check_payload_identity(raw_path: Path, canonical_path: Path, witness: dict) -> tuple[dict, dict]:
    for path, field in ((raw_path, "raw_payload_sha256"), (canonical_path, "canonical_payload_sha256")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != witness[field]:
            raise ValueError(f"retained witness identity changed: {path.name}: {digest} != {witness[field]}")
    raw, canonical = BASE.read_payload(raw_path), BASE.read_payload(canonical_path)
    if (raw["mode_dim"], raw["t0"], raw["t1"]) != (witness["dimension"], witness["word_t0"], witness["word_t1"]):
        raise ValueError("retained source mode or endpoints changed")
    for a, b in zip(raw["events"], canonical["events"], strict=True):
        if (a["type"], a["time"]) != (b["type"], b["time"]):
            raise ValueError("canonicalization detached the event history")
        if a["type"] == BASE.EV_S and not np.array_equal(a["R"], b["R"]):
            raise ValueError("actual applied R_S was changed")
    return raw, canonical


def audit_mode(raw_path: Path, canonical_path: Path, witness: dict,
               domain: dict, source_input: Path, trace=None) -> dict:
    raw, payload = check_payload_identity(raw_path, canonical_path, witness)
    radius = domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"]
    linear = BASE.build_reset_normalized_linear_path(payload, radius)
    # Do not use the freshly computed maximizing direction or search new scales.
    linear["direction"] = np.asarray(witness["direction"], dtype=float)
    if linear["counts"] != witness["event_counts"]:
        raise ValueError("retained complete-word operation counts changed")
    kernels = FAST.prepare_kernels(payload, linear)
    c = BIAS.build()["BIAS0"]
    pred_checks = []
    cov_differences = []
    for ev, raw_ev, lp in zip(payload["events"], raw["events"], linear["path"], strict=True):
        cov_differences.append(float(np.linalg.norm(lp["Pafter"] - raw_ev["Pafter_shipping"][:payload["mode_dim"], :payload["mode_dim"]]) / max(1.0, np.linalg.norm(raw_ev["Pafter_shipping"]))))
        if payload["mode_dim"] == 21 and ev["type"] == BASE.EV_PRED:
            expected = float(np.float32(math.exp(-ev["h"] / ev["tau_ba"])))
            pred_checks.append(ev["tau_ba"] == c["filter_tau_s"] and np.array_equal(ev["linear_shipping"][18:21, 18:21], expected * np.eye(3)))
    if pred_checks and not all(pred_checks):
        raise ValueError("retained bias prediction no longer matches the shipping GM factor")
    cases = []
    for old in witness["cases"]:
        case = audit_case(payload, linear, kernels, domain, float(old["scale"]), source_input, trace)
        case["retained_report_rho_endpoint"] = old["rho_endpoint"]
        case["rho_reproduction_difference"] = case["rho_endpoint"] - old["rho_endpoint"]
        if abs(case["rho_reproduction_difference"]) > 1e-8:
            raise ValueError("fixed retained state no longer reproduces its point endpoint energy")
        if not case["bias_recurrence_point_check_pass"]:
            raise ValueError("retained finite evaluator disagrees with BIAS1 corrected-error recurrence")
        cases.append(case)
    passing = [x for x in cases if x["conditional_point_checks_pass"]]
    return {
        "dimension": payload["mode_dim"], "word_t0": payload["t0"], "word_t1": payload["t1"],
        "event_counts": linear["counts"], "raw_payload_sha256": witness["raw_payload_sha256"],
        "canonical_payload_sha256": witness["canonical_payload_sha256"],
        "retained_payload_identity_verified": True,
        "all_actual_RS_values_copied_bit_for_bit": True,
        "shipping_float32_GM_factors_verified": all(pred_checks),
        "GM_factor_rounding_is_numerical_not_a_new_physical_premise": True,
        "source_nonzero_attitude_reset_count": sum(bool(np.linalg.norm(e["dtheta"]) > 0) for e in raw["events"]),
        "max_reset_normalized_vs_captured_covariance_relative_difference": max(cov_differences),
        "source_evidence": ATTACH.source_evidence_inventory(raw),
        "captured_reset_storage_transport": ATTACH.audit_captured_transport(raw, linear),
        "cases": cases,
        "cases_passing_conditional_point_checks": len(passing),
        "expanding_cases_passing_conditional_point_checks": [x["scale"] for x in passing if x["rho_endpoint"] > 1.0],
        "worst_rho_passing_conditional_point_checks": max((x["rho_endpoint"] for x in passing), default=None),
        "minimum_BIAS2_sufficient_point_ratio_on_checked_cases": min((x["BIAS2"]["mu_sufficient_point"] for x in passing), default=None),
        "this_sample_minimum_is_not_a_uniform_separation_constant": True,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload-prefix", type=Path, required=True)
    ap.add_argument("--witnesses", type=Path, default=DEFAULT_WITNESSES)
    ap.add_argument("--domain", type=Path, default=BIAS.DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--unpack-retained", action="store_true")
    args = ap.parse_args()
    witness = json.loads(args.witnesses.read_text())
    if args.unpack_retained:
        unpack_retained_payloads(args.witnesses.with_name("ou3_p4_retained_payloads.tar.xz"), args.payload_prefix, witness)
    domain = json.loads(args.domain.read_text())
    if witness["true_bias_root"] != [0.0] * 3 or witness["deterministic_bias_mismatch"] != [0.0] * 3:
        raise ValueError("this audit must retain the archived homogeneous zero-bias source")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    trace_path = args.output.with_suffix(".prefixes.csv.gz")
    modes = {}
    with gzip.open(trace_path, "wt", newline="") as stream:
        trace = csv.DictWriter(stream, fieldnames=TRACE_FIELDS)
        trace.writeheader()
        for mode in ("H18", "A21"):
            modes[mode] = audit_mode(
                Path(f"{args.payload_prefix}_{mode}.bin"),
                Path(f"{args.payload_prefix}_canonical_{mode}.bin"),
                witness["modes"][mode], domain, Path(witness["source_input_name"]), trace,
            )
    report = {
        "qualification": "NON_PROMOTING_RETAINED_BIAS_WITNESS_ADMISSIBILITY_AND_BIAS2_POINT_AUDIT",
        "conditional_mathematical_experiment": True,
        "hardware_qualification_required_to_run_this_experiment": False,
        "hardware_qualification_closed_here": False,
        "new_source_or_direction_or_scale_search": False,
        "filter_or_domain_changed": False,
        "BIAS2_evaluated_on_actual_corrected_error_histories": True,
        "BIAS2_uniform_sector_available": False,
        "source_membership_status": "UNRESOLVED: archived payload lacks a complete BRMM joint source/bias realization witness",
        "finite_reset_storage_attachment_status": "REJECTED_AS_GAUGE_ATTACHMENT: deleting captured resets without transporting F/Q/H does not attach this finite word; exact finite pullback and signed reset identities remain valid",
        "critic": "An admitted rho>1 point cannot satisfy a negative master after adding nonnegative multiples of sectors valid at that point, even when its BIAS2 separation ratio is positive.",
        "canonical_P4_falsified_here": False, "P4_promoted": False, "P5_may_start": False,
        "modes": modes, "prefix_trace": trace_path.name,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    for mode, data in modes.items():
        print("BIAS_WITNESS_AUDIT", mode, json.dumps({k: data[k] for k in (
            "cases_passing_conditional_point_checks", "expanding_cases_passing_conditional_point_checks",
            "worst_rho_passing_conditional_point_checks", "minimum_BIAS2_sufficient_point_ratio_on_checked_cases",
            "max_reset_normalized_vs_captured_covariance_relative_difference",
        )}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
