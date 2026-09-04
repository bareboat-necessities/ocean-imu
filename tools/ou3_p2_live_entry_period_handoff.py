#!/usr/bin/env python3
"""Bind WavePeriodEstimator startup takeover to the Normal-Live P2 source language.

The deployed tuner now receives a usable WavePeriodEstimator estimate before the
filter may enter Live.  This certificate makes the startup-to-Live handoff
explicit for the P2/P3 source language:

* before startup qualification, tuner input may be the fixed 0.2 Hz prior;
* WavePeriodEstimator::hasUsablePeriod() is a one-way latch from the same
  canonical period state;
* TunerReady/Live require that latch, but not strict isReady();
* the tuner consumes the previous sample's period state, so the sample that
  first latches usable may still have staged a candidate using the prior;
* consequently the initial Normal-Live committed/staged pair must retain
  pre-Live prior influence, while every *future* tuner frequency selection in
  Normal Live comes from the estimator path.

The existing P2 tuner graph is deliberately broader than this handoff: it keeps
the full clamped frequency interval and does not prune the initial pair.  That
is source-safe.  A numeric value equal to 0.2 Hz remains in the graph because
the estimator may itself produce that value; this must not be confused with
re-enabling the fixed-prior selector branch after Live entry.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_source_path_reachability as PATH
import ou3_sea3_wave_period_frontend as FRONT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
SCHEMA = 1
QUALIFICATION = "OU3_P2_LIVE_ENTRY_WAVE_PERIOD_HANDOFF"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    if path != DEFAULT_DOMAIN.resolve():
        raise ValueError(
            "live-entry handoff is bound to the repository operating-domain artifact"
        )

    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("live-entry handoff must not be trajectory fitted")

    front = FRONT.build(REPO)
    ff = FRONT.validate(front)
    startup = front.get("startup_source_language", {})
    parity = front.get("source_parity", {})
    domain_parity = front.get("operating_domain_parity", {})

    estimator_text = ESTIMATOR.read_text(encoding="utf-8")
    compact_estimator = " ".join(estimator_text.split())
    usable_implies_finite_positive_period = all(
        " ".join(marker.split()) in compact_estimator
        for marker in (
            "const float period = getPeriodSec();",
            "if (!(std::isfinite(period) && period > 0.0f)) return;",
            "usable_period_ = true;",
        )
    )

    constants = PATH._constants()
    f_min = float(constants["min_freq"])
    f_max = float(constants["max_freq"])
    prior = float(domain["startup"]["tune_frequency_prior_hz"])
    screened = list(map(float, front["declared_inputs"]["screened_tuning_frequency_hz"]))
    screened_contains_graph_clamp = (
        len(screened) == 2
        and screened[0] <= f_min <= f_max <= screened[1]
    )
    prior_inside_graph_numeric_range = f_min <= prior <= f_max

    live_requires_usable = startup.get(
        "live_entry_requires_wave_period_estimator_usable"
    ) is True
    takeover_latched = startup.get(
        "wave_period_startup_takeover_is_one_way_latched"
    ) is True
    selector_requires_usable = parity.get(
        "selector_requires_startup_usable_period"
    ) is True
    selector_falls_back_to_prior = parity.get(
        "selector_falls_back_to_fixed_prior"
    ) is True
    strict_ready_not_required = (
        startup.get("tuner_ready_requires_wave_period_estimator_ready") is False
        and startup.get("wave_period_takeover_waits_for_isReady") is False
    )

    fixed_prior_selector_reachable_in_normal_live = not (
        live_requires_usable
        and takeover_latched
        and selector_requires_usable
        and selector_falls_back_to_prior
        and usable_implies_finite_positive_period
    )

    failures = list(ff)
    if not all(domain_parity.values()):
        failures.append("operating-domain startup fields do not match source defaults")
    if startup.get("tuner_consumes_previous_sample_wave_period_state") is not True:
        failures.append("tuner/period one-sample ordering certificate disappeared")
    if startup.get(
        "first_newly_usable_wave_period_can_affect_tuner_no_earlier_than_next_valid_sample"
    ) is not True:
        failures.append("first-usable one-sample takeover edge disappeared")
    if fixed_prior_selector_reachable_in_normal_live:
        failures.append("fixed-prior selector branch remains reachable in Normal Live")
    if not screened_contains_graph_clamp:
        failures.append("SEA0 screened frequency interval no longer contains P2 clamp")
    if not prior_inside_graph_numeric_range:
        failures.append("fixed startup prior left the P2 numeric frequency range")
    if not strict_ready_not_required:
        failures.append("startup unexpectedly requires strict WavePeriodEstimator readiness")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "wave_period_frontend_schema": front.get("schema_version"),
        "wave_period_frontend_validation_pass": not ff,
        "live_entry_requires_wave_period_estimator_usable": live_requires_usable,
        "live_entry_requires_wave_period_estimator_strict_ready": False,
        "wave_period_usable_takeover_one_way_latched": takeover_latched,
        "usable_period_implies_finite_positive_period": usable_implies_finite_positive_period,
        "tuner_consumes_previous_sample_wave_period_state": startup.get(
            "tuner_consumes_previous_sample_wave_period_state"
        ) is True,
        "first_usable_takeover_has_one_valid_sample_delay": startup.get(
            "first_newly_usable_wave_period_can_affect_tuner_no_earlier_than_next_valid_sample"
        ) is True,
        "normal_live_fixed_prior_selector_branch_reachable":
            fixed_prior_selector_reachable_in_normal_live,
        "normal_live_future_frequency_source":
            "WavePeriodEstimator canonical startup-usable state, then source clamp",
        "p2_frequency_clamp_hz": [f_min, f_max],
        "fixed_prior_frequency_hz": prior,
        "fixed_prior_numeric_value_inside_p2_frequency_range":
            prior_inside_graph_numeric_range,
        "numeric_prior_value_does_not_imply_prior_selector_branch": True,
        "pre_live_prior_influence_retained_in_initial_live_pair": True,
        "initial_live_committed_staged_pair_pruned_by_this_certificate": False,
        "initial_live_pair_reason": (
            "the usable latch may occur after the current sample tuner update; "
            "a candidate staged before takeover may therefore remain in the "
            "committed/staged state at the Live boundary"
        ),
        "existing_p2_full_frequency_range_is_conservative_for_future_live_staging":
            screened_contains_graph_clamp,
        "strict_isReady_remains_diagnostic_only_for_live_entry": strict_ready_not_required,
        "P2_LIVE_ENTRY_WAVE_PERIOD_HANDOFF_CERTIFICATE":
            "PASS" if not failures else "FAIL",
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "next_obligation": (
            "retain the unrestricted source-safe initial Live pair, but do not "
            "reintroduce the fixed-prior selector branch inside Normal-Live words; "
            "future SEA3 pruning may use finite EW/log-period dynamics once certified"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    for key in (
        "source_only",
        "wave_period_frontend_validation_pass",
        "live_entry_requires_wave_period_estimator_usable",
        "wave_period_usable_takeover_one_way_latched",
        "usable_period_implies_finite_positive_period",
        "tuner_consumes_previous_sample_wave_period_state",
        "first_usable_takeover_has_one_valid_sample_delay",
        "numeric_prior_value_does_not_imply_prior_selector_branch",
        "pre_live_prior_influence_retained_in_initial_live_pair",
        "existing_p2_full_frequency_range_is_conservative_for_future_live_staging",
        "strict_isReady_remains_diagnostic_only_for_live_entry",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "live_entry_requires_wave_period_estimator_strict_ready",
        "normal_live_fixed_prior_selector_branch_reachable",
        "initial_live_committed_staged_pair_pruned_by_this_certificate",
        "P3_PROMOTED",
        "P4_PROMOTED",
        "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    clamp = d.get("p2_frequency_clamp_hz", [])
    if (
        not isinstance(clamp, list)
        or len(clamp) != 2
        or not all(isinstance(x, (int, float)) and math.isfinite(float(x)) for x in clamp)
        or not (0.0 < float(clamp[0]) < float(clamp[1]))
    ):
        f.append("invalid P2 frequency clamp")
    if d.get("P2_LIVE_ENTRY_WAVE_PERIOD_HANDOFF_CERTIFICATE") != "PASS":
        f.append("live-entry WavePeriodEstimator handoff did not pass")
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
    print(json.dumps({
        "status": d["P2_LIVE_ENTRY_WAVE_PERIOD_HANDOFF_CERTIFICATE"],
        "live_requires_usable": d["live_entry_requires_wave_period_estimator_usable"],
        "prior_selector_reachable_in_live": d[
            "normal_live_fixed_prior_selector_branch_reachable"
        ],
        "prelive_prior_influence_retained": d[
            "pre_live_prior_influence_retained_in_initial_live_pair"
        ],
        "frequency_clamp_hz": d["p2_frequency_clamp_hz"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
