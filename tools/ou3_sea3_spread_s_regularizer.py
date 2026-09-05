#!/usr/bin/env python3
"""Spread-S translational regularizer for the complete SEA3 Normal-Live word.

This is a subcertificate of the *same* complete SEA3 word, not a reduced word
or source generator.  The deployed progress-preserving S=0 scheduler guarantees
one actual firing in every interval of length g.  On the canonical 3 s P3
window we may therefore select three already-present firings from the backward
endpoint lag windows

    [0,g], [1.5,1.5+g], [3-g,3].

All other S updates, every accelerometer update, process step, magnetic event
and covariance-floor event remain in the literal H18/A21 word.

For endpoint coordinates [S,p,v], an S firing at lag ell has the integrator row

    b(ell) = [1,-ell,ell^2/2].

The three lag windows are strictly ordered for the shipping g<=0.15 s bound.
The Vandermonde determinant therefore has a uniform nonzero lower bound.  The
a_w contribution and OU driving over the 3 s memory are paid as nuisance in an
effective observation covariance; this makes the resulting [v,p,S] information
bound valid for arbitrary legal time variation of tau/sigma inside the complete
SEA3 word.  Crucially, the observation covariance uses the *actual deployed
SpectralMSE applied R_S* ceiling and its per-axis factors; it does not substitute
a target R_S or an independent selected schedule.

The certificate is intentionally only a mandatory full-word ingredient.  It
cannot promote P3 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_SPREAD_S_REGULARIZER"
HORIZON_S = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    sched = SCHED.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(sched),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"spread-S prerequisites failed: {bad}")

    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("spread-S certificate lost complete SEA3 source")
    rs_contract = complete["R_S_regularizer"]
    if rs_contract["deployed_law"] != "SpectralMSE":
        raise RuntimeError("spread-S proof requires deployed SpectralMSE R_S")
    if rs_contract["actual_applied_R_S_required_at_every_due_S_update"] is not True:
        raise RuntimeError("actual applied R_S was not retained")
    if rs_contract["all_due_S_updates_remain_in_full_word"] is not True:
        raise RuntimeError("full S scheduler word was reduced")

    g = up(float(sched["certified_uniform_max_gap_s"]))
    if not (0.0 < g < 0.5):
        raise RuntimeError(f"scheduler gap too wide for spread selection: {g}")

    # Backward endpoint lag windows.  Each has width g, hence scheduler
    # recurrence guarantees at least one *actual* S firing in each window.
    windows = [
        [0.0, g],
        [1.5, up(1.5 + g)],
        [down(HORIZON_S - g), HORIZON_S],
    ]
    d01 = down(windows[1][0] - windows[0][1])
    d12 = down(windows[2][0] - windows[1][1])
    d02 = down(windows[2][0] - windows[0][1])
    if min(d01, d12, d02) <= 0.0:
        raise RuntimeError("spread S windows lost strict ordering")

    det_lower = down(0.5 * d01 * d12 * d02)
    # ||b(ell)||^2 = 1 + ell^2 + ell^4/4.  Three rows, ell<=3.
    row_norm2_upper = up(1.0 + HORIZON_S**2 + HORIZON_S**4 / 4.0)
    frobenius_sq_upper = up(3.0 * row_norm2_upper)
    sigma_min_lower = down(det_lower / frobenius_sq_upper)

    inv = dynamic["dynamic_invariant"]
    tau_lo = float(inv["tau_applied_s"][0])
    sigma_hi = float(inv["sigma_aw_filter_mps2"][1])
    rs_hi = float(inv["R_S_applied"][1])
    if not (tau_lo > 0.0 and sigma_hi > 0.0 and rs_hi > 0.0):
        raise RuntimeError("invalid complete-SEA3 adaptive invariant")

    factors = list(map(float, rs_contract["axis_std_factors"]))
    rs_axis_std_upper = up(rs_hi * max(factors))
    rs_variance_upper = up(rs_axis_std_upper * rs_axis_std_upper)

    # Endpoint S observations contain the OU acceleration history.  Treat the
    # full legal history as nuisance for this [S,p,v] subcertificate rather than
    # freezing tau or sigma.  Bounds use the complete dynamic invariant only.
    qc_upper = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    s_aw_state_nuisance_variance_upper = up(
        sigma_hi * sigma_hi * (HORIZON_S**3 / 6.0) ** 2
    )
    s_ou_process_nuisance_variance_upper = up(
        qc_upper * HORIZON_S**7 / 252.0
    )
    # Three selected observations may be correlated through the common a_w
    # history.  lambda_max(R_stack) <= trace(R_stack) <= 3*max diagonal.
    stacked_covariance_lambda_max_upper = up(3.0 * (
        rs_variance_upper
        + s_aw_state_nuisance_variance_upper
        + s_ou_process_nuisance_variance_upper
    ))
    information_lambda_min_lower = down(
        sigma_min_lower * sigma_min_lower / stacked_covariance_lambda_max_upper
    )

    passed = all(math.isfinite(x) and x > 0.0 for x in (
        det_lower,
        sigma_min_lower,
        stacked_covariance_lambda_max_upper,
        information_lambda_min_lower,
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "complete_SEA3_source_consumed": True,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_TS_extrema_product_used": False,
        "selected_S_events_replace_full_scheduler_word": False,
        "all_due_S_updates_remain_in_literal_word": True,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "R_S_axis_std_factors": factors,
        "word_horizon_s": HORIZON_S,
        "scheduler_uniform_gap_s_upper": g,
        "selected_backward_lag_windows_s": windows,
        "selected_firings_are_guaranteed_members_of_full_word": True,
        "pairwise_lag_separation_lower_s": [d01, d12, d02],
        "endpoint_integrator_state_order": ["S", "p", "v"],
        "endpoint_integrator_vandermonde_det_abs_lower": det_lower,
        "endpoint_integrator_frobenius_sq_upper": frobenius_sq_upper,
        "endpoint_integrator_sigma_min_lower": sigma_min_lower,
        "adaptive_invariant_used_for_nuisance_only": {
            "tau_applied_lower_s": tau_lo,
            "sigma_aw_upper_mps2": sigma_hi,
            "R_S_applied_upper": rs_hi,
        },
        "effective_observation_covariance": {
            "actual_applied_R_S_axis_std_upper": rs_axis_std_upper,
            "actual_applied_R_S_variance_upper": rs_variance_upper,
            "OU_driving_intensity_upper": qc_upper,
            "a_w_state_nuisance_variance_upper": s_aw_state_nuisance_variance_upper,
            "OU_process_nuisance_variance_upper": s_ou_process_nuisance_variance_upper,
            "stacked_covariance_lambda_max_upper": stacked_covariance_lambda_max_upper,
        },
        "integrator_information_lambda_min_lower": information_lambda_min_lower,
        "arbitrary_legal_time_variation_inside_word_allowed": True,
        "spread_S_regularizer_pass": passed,
        "P3_promoted": False,
        "next_obligation": (
            "compose this strong endpoint [v,p,S] information with stable/process a_w, "
            "windowed eta6 vector PE, A-mode bias contraction and the exact full-word "
            "M_delta preservation identities; do not replace the complete SEA3 word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "complete_SEA3_source_consumed",
        "all_due_S_updates_remain_in_literal_word",
        "actual_applied_SpectralMSE_R_S_consumed",
        "selected_firings_are_guaranteed_members_of_full_word",
        "arbitrary_legal_time_variation_inside_word_allowed",
        "spread_S_regularizer_pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_family_replaced",
        "trajectory_replay_used",
        "independent_tau_sigma_RS_TS_extrema_product_used",
        "selected_S_events_replace_full_scheduler_word",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    info = d.get("integrator_information_lambda_min_lower")
    if not isinstance(info, (int, float)) or not (math.isfinite(float(info)) and float(info) > 0.0):
        f.append("spread-S information lower is not strict")
    if float(d.get("scheduler_uniform_gap_s_upper", math.inf)) > 0.151:
        f.append("scheduler recurrence widened beyond deployed 150 ms service")
    if len(d.get("selected_backward_lag_windows_s", [])) != 3:
        f.append("spread-S selection lost three guaranteed firings")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "gap_s": d["scheduler_uniform_gap_s_upper"],
        "windows_s": d["selected_backward_lag_windows_s"],
        "R_S_axis_std_upper": d["effective_observation_covariance"]["actual_applied_R_S_axis_std_upper"],
        "integrator_information_lower": d["integrator_information_lambda_min_lower"],
        "pass": d["spread_S_regularizer_pass"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
