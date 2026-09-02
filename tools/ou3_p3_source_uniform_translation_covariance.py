#!/usr/bin/env python3
"""Source-uniform directional translation covariance upper for OU-III P3.

The existing finite-memory estimator in :mod:`ou3_source_reachable_matrix_p3`
was historically evaluated one source cell at a time.  Its inequalities do not
require the tuner parameters to stay in one cell, however: every term is built
from pointwise source extrema over the observation window.

Evaluating it on the complete invariant source box therefore gives a single
Loewner-safe diagonal dominator for [v,p,S,a_w] under arbitrary normal-Live
variation of tau, sigma_aw and R_S:

* pseudo-update gaps use the largest cadence admitted by the full tau box;
* OU nuisance/process terms use sigma_max and tau_min;
* selected S observations use the largest deployed R_S standard deviation;
* endpoint propagation uses the resulting source-uniform word horizon.

No Cartesian *sequence* is asserted: the extrema are used only as one-sided
bounds in a finite-memory estimator.  Allowing them simultaneously can enlarge
the covariance upper but cannot exclude a source trajectory.  This producer is
therefore source complete but intentionally conservative.

It supplies only the covariance upper.  P3 still needs a source-uniform
post-measurement process lower before theorem promotion.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("source-uniform covariance upper must not be trajectory fitted")
    live = domain["normal_live"]
    sched = BASE.source_schedule()
    tau = Interval(*map(float, sched["tau_applied_invariant_s"]))
    sigma = Interval(*map(float, sched["sigma_aw_applied_safety"]))
    rs = Interval(*map(float, sched["R_S_applied_invariant"]))
    Tpe = BASE.pos(live["vector_pe_recurrence_window_s"], "PE recurrence")
    upper, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)

    failures = []
    if len(upper) != 4 or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in upper):
        failures.append("translation covariance upper is not four finite positive entries")
    if not float(timing["word_horizon_s_lower"]) > 0.0:
        failures.append("word horizon lower is not positive")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_SOURCE_UNIFORM_TIME_VARYING_TRANSLATION_COVARIANCE_UPPER",
        "source_generated_not_trajectory_fit": True,
        "time_varying_source_parameters_covered_by_pointwise_extrema": True,
        "finite_memory_estimator_route": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "tau_applied_invariant_s": tau.as_list(),
        "sigma_aw_applied_safety_mps2": sigma.as_list(),
        "R_S_applied_invariant_std": rs.as_list(),
        "state_order": ["v", "p", "S", "a_w"],
        "Sigma_translation_diagonal_upper": [float(x) for x in upper],
        "timing": timing,
        "P3_PROCESS_LOWER_ESTABLISHED_HERE": False,
        "P3_PROMOTED": False,
        "next_obligation": (
            "compare the source-uniform time-varying process floor, after all interleaved measurements, "
            "against this directional covariance upper"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_SOURCE_UNIFORM_TIME_VARYING_TRANSLATION_COVARIANCE_UPPER":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "time_varying_source_parameters_covered_by_pointwise_extrema",
        "finite_memory_estimator_route",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed",
        "P3_PROCESS_LOWER_ESTABLISHED_HERE", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("state_order") != ["v", "p", "S", "a_w"]:
        f.append("translation state order changed")
    upper = d.get("Sigma_translation_diagonal_upper", [])
    if len(upper) != 4 or any(not (isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) > 0.0) for x in upper):
        f.append("invalid translation covariance upper")
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
        "Sigma_translation_diagonal_upper": d["Sigma_translation_diagonal_upper"],
        "word_horizon_s_lower": d["timing"]["word_horizon_s_lower"],
        "word_horizon_s_upper": d["timing"]["word_horizon_s_upper"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
