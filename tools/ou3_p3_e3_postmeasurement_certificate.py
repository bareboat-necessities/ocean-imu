#!/usr/bin/env python3
"""e3-sharpened source-uniform post-measurement P3 certificate.

This certificate keeps the source-complete covariance, clock-phase decay,
selected process modes, and one-shot lifted measurement conditioning of the
retained P3 LTV route.  It changes only the spectral conversion of the proven
4-D translation Gramian determinant: :mod:`ou3_p3_gramian_e3` replaces the
very loose ``det/trace^3`` step by the rigorous ``det/e3`` bound.

The producer is intentionally fail-closed.  It records the gain from that
mathematical sharpening but promotes canonical P3 only if both H and A remain
above the unchanged 1e-18 usefulness threshold after the same measurement
attenuation used by the retained route.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_gramian_e3 as E3
import ou3_p3_ltv_postmeasurement_certificate as POST
import ou3_p3_ltv_translation_ucc_probe as LTV
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _mode_certificate(pre: dict, mode: str, blocks: dict, domain: dict) -> dict:
    upper = [
        float(x)
        for x in pre["source_complete_translation"]["covariance_upper"]
        ["Sigma_translation_diagonal_upper"]
    ]
    rows = []
    for endpoint in pre["source_complete_translation"]["endpoint_rows"]:
        candidates = []
        for probe in endpoint["candidates"]:
            sharp = E3.sharpen_probe(probe, upper)
            row = POST._info_attenuation(sharp, mode, blocks, domain)
            row["endpoint_tau_index"] = int(endpoint["endpoint_tau_index"])
            row["clock_phase_decay_exponent_upper"] = float(
                probe["decay_exponent_upper"]
            )
            row["trace3_translation_floor_baseline"] = float(
                sharp["trace3_relative_process_floor_lower_baseline"]
            )
            row["e3_translation_floor_lower"] = float(
                sharp["relative_process_floor_lower"]
            )
            row["e3_over_trace3_translation_improvement"] = float(
                sharp["improvement_over_trace3"]
            )
            row["Sigma_normalized_gramian_e3_upper"] = float(
                sharp["Sigma_normalized_gramian_e3_upper"]
            )
            candidates.append(row)
        best = max(
            candidates,
            key=lambda x: x["relative_Riccati_injection_margin_lower"],
        )
        rows.append(
            {
                "endpoint_tau_index": int(endpoint["endpoint_tau_index"]),
                "best": best,
                "candidates": candidates,
            }
        )
    worst = min(
        rows,
        key=lambda x: x["best"]["relative_Riccati_injection_margin_lower"],
    )
    delta = float(worst["best"]["relative_Riccati_injection_margin_lower"])
    return {
        "mode": mode,
        "endpoint_tau_cells_scanned": len(rows),
        "endpoint_rows": rows,
        "worst_endpoint": worst,
        "relative_Riccati_injection_margin_lower": delta,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "useful_margin_established": delta >= BASE.MIN_USEFUL_DELTA,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("e3 P3 certificate must not be trajectory fitted")
    POST._source_contract(domain)

    pre = LTV.build(path, source_node_indices=())
    pf = LTV.validate(pre)
    if pf:
        raise RuntimeError(f"retained LTV source producer failed: {pf}")
    blocks = POST._source_uniform_endpoint_blocks(domain)
    modes = {m: _mode_certificate(pre, m, blocks, domain) for m in ("H", "A")}
    worst = min(
        float(modes[m]["relative_Riccati_injection_margin_lower"])
        for m in ("H", "A")
    )
    passed = all(modes[m]["useful_margin_established"] for m in ("H", "A"))

    old = POST.build(path)
    of = POST.validate(old)
    if of:
        raise RuntimeError(f"retained postmeasurement baseline failed validation: {of}")
    old_worst = float(old["worst_H_A_relative_Riccati_injection_margin_lower"])
    improvement = math.inf if old_worst <= 0.0 else worst / old_worst

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_SOURCE_UNIFORM_E3_LTV_POSTMEASUREMENT_CERTIFICATE",
        "source_generated_not_trajectory_fit": True,
        "linear_P3_only": True,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "retained_LTV_determinant_consumed": True,
        "trace_cubed_spectral_conversion_used": False,
        "e3_hadamard_spectral_conversion_used": True,
        "numerical_eigendecomposition_used": False,
        "same_lifted_measurement_attenuation_as_retained_route": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "retained_trace3_worst_H_A_margin": old_worst,
        "e3_worst_H_A_margin": worst,
        "improvement_over_retained_trace3_margin": improvement,
        "modes": modes,
        "P3_LINEAR_CERTIFICATE_ESTABLISHED": passed,
        "P3_PROMOTED": passed,
        "next_obligation": (
            "feed the promoted e3 P3 metric into the finite-angle P4 complete-word route"
            if passed
            else "retain source/path correlation in covariance and measurement bounds; the e3 spectral loss is no longer the dominant bottleneck"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_SOURCE_UNIFORM_E3_LTV_POSTMEASUREMENT_CERTIFICATE":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "linear_P3_only",
        "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch",
        "retained_LTV_determinant_consumed",
        "e3_hadamard_spectral_conversion_used",
        "same_lifted_measurement_attenuation_as_retained_route",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trace_cubed_spectral_conversion_used",
        "numerical_eigendecomposition_used",
        "trajectory_replay_used",
        "filter_changed",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    old = d.get("retained_trace3_worst_H_A_margin")
    new = d.get("e3_worst_H_A_margin")
    if not isinstance(old, (int, float)) or not (math.isfinite(float(old)) and float(old) > 0.0):
        f.append("retained trace3 margin is not strict")
    if not isinstance(new, (int, float)) or not (math.isfinite(float(new)) and float(new) > 0.0):
        f.append("e3 margin is not strict")
    if isinstance(old, (int, float)) and isinstance(new, (int, float)) and float(new) <= float(old):
        f.append("e3 spectral bound did not improve the retained trace3 margin")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        x = m.get("relative_Riccati_injection_margin_lower")
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"{mode}: e3 P3 margin is not strict")
        if int(m.get("endpoint_tau_cells_scanned", 0)) != 10:
            f.append(f"{mode}: did not scan ten endpoint tau cells")
    passed = all(
        float(d["modes"][m]["relative_Riccati_injection_margin_lower"])
        >= BASE.MIN_USEFUL_DELTA
        for m in ("H", "A")
    ) if all(m in d.get("modes", {}) for m in ("H", "A")) else False
    if d.get("P3_PROMOTED") is not passed:
        f.append("P3 promotion flag does not match unchanged usefulness gate")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(
        json.dumps(
            {
                "retained_trace3_worst": d["retained_trace3_worst_H_A_margin"],
                "e3_worst": d["e3_worst_H_A_margin"],
                "improvement": d["improvement_over_retained_trace3_margin"],
                "P3_PROMOTED": d["P3_PROMOTED"],
                "validation_failures": vf,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
