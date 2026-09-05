#!/usr/bin/env python3
"""Canonical hard finite-window SEA3 provider contract.

This module owns the *only* artifact shape that may materialize the 3 s source
family consumed by canonical P3.  It does not create a second source.  The
provider must certify one compact, phase-continuous SEA3 realization family

    zeta_k = (x^s_k, lambda_k, z^t_k, q_k)

for all 601 samples.  The same x^s/lambda history must generate the physical
translation/rotation, the private measurement-only front end, the tuner and
scheduler, and hence the complete H18/A21 Riccati word.

The SEA3 theorem permits either an oscillator/shaping realization or an
equivalent hard finite-window dynamic constraint.  Both are represented by the
same transition-witness contract below.  A replay, a Gaussian good event,
spectral moments alone, arbitrary per-sample boxes, a fixed-lambda word, a
finite RAO grid, or an independently selected tuner schedule cannot satisfy
this provider.

SEA3 compactness is already a theorem-domain fact.  The open obligation here is
strictly machine execution: the repository does not yet encode a validated
finite-window realization of the theorem's compact transition relation and
shaping state.  In particular, the current continuum RAO producer is a set/
moment theorem, not an x^s trajectory generator.  This module therefore keeps
three separate executable-provider obligations visible without reinterpreting
any of them as a compactness failure:

* machine-readable validated realization of R_lambda over the 601 samples;
* hard shaping-state/excitation (or equivalent hard-window IQC) propagation;
* joint translational/rotational response from that same shaping history.

``PROVIDER_IMPLEMENTATION_CLOSED`` is code-owned and remains false until those
obligations are implemented and validated.  An input JSON cannot self-assert
its way into P3.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_directional_response_family as RESPONSE
import ou3_sea3_physical_admissibility as PHYSICAL

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_HARD_FINITE_WINDOW_REALIZATION_V2"
CANONICAL_SOURCE = "COMPLETE_SEA3_NORMAL_LIVE_WORD"
HORIZON_S = 3.0
DT_S = 0.005
SAMPLES = 601

# These gates describe executable proof machinery, not SEA3 set compactness.
MACHINE_READABLE_R_LAMBDA_CLOSED = False
HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED = False
JOINT_TRANSLATIONAL_ROTATIONAL_SHAPING_CLOSED = False
PROVIDER_IMPLEMENTATION_CLOSED = (
    MACHINE_READABLE_R_LAMBDA_CLOSED
    and HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED
    and JOINT_TRANSLATIONAL_ROTATIONAL_SHAPING_CLOSED
)

_FORBIDDEN_TRUE_FLAGS = (
    "trajectory_replay_used",
    "gaussian_good_event_used",
    "spectral_moment_only_source_used",
    "arbitrary_bounded_input_source_used",
    "fixed_lambda_word_used",
    "independent_axis_boxes_used",
    "independent_SEA_RAO_cartesian_product_used",
    "finite_RAO_grid_used",
    "independent_tuner_schedule_used",
    "retired_P2_graph_used",
    "selected_four_S_word_used",
)

_REQUIRED_TRUE_FLAGS = (
    "SEA3_parameter_domain_compact",
    "compact_transition_relation_is_theorem_domain",
    "phase_continuous",
    "same_xs_lambda_history_for_all_channels",
    "same_realization_drives_translation_rotation_frontend_tuner_geometry",
    "all_valid_accelerometer_samples_retained",
    "all_due_S_updates_retained",
    "actual_applied_per_axis_RS_retained",
    "asynchronous_vector_PE_events_retained",
    "covariance_floor_events_retained",
)


def _continuity_failures(samples: Any) -> list[str]:
    f: list[str] = []
    if not isinstance(samples, list) or len(samples) != SAMPLES:
        return [f"window must contain exactly {SAMPLES} source transitions"]
    previous: dict[str, Any] | None = None
    source_id: str | None = None
    for k, sample in enumerate(samples):
        if not isinstance(sample, dict):
            f.append(f"sample {k} is not an object")
            continue
        if sample.get("k") != k:
            f.append(f"sample {k} index mismatch")
        sid = sample.get("source_identity")
        if not isinstance(sid, str) or not sid:
            f.append(f"sample {k} missing source_identity")
        elif source_id is None:
            source_id = sid
        elif sid != source_id:
            f.append(f"sample {k} changed source_identity")
        for key in (
            "xs_in_id", "xs_out_id", "lambda_in_id", "lambda_out_id",
            "source_transition_witness_id", "joint_response_witness_id",
        ):
            if not isinstance(sample.get(key), str) or not sample[key]:
                f.append(f"sample {k} missing {key}")
        if previous is not None:
            if previous.get("xs_out_id") != sample.get("xs_in_id"):
                f.append(f"sample {k} broke x^s phase continuity")
            if previous.get("lambda_out_id") != sample.get("lambda_in_id"):
                f.append(f"sample {k} broke lambda transition continuity")
        for key in ("joint_physical_output", "source_events"):
            if not isinstance(sample.get(key), dict):
                f.append(f"sample {k} missing {key}")
        previous = sample
    return f


def validate_candidate_structure(d: dict[str, Any]) -> list[str]:
    """Validate artifact shape and anti-shortcut semantics, but not SEA0 math."""
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("wrong canonical source")
    if not math.isclose(float(d.get("window_horizon_s", math.nan)), HORIZON_S, rel_tol=0.0, abs_tol=1e-12):
        f.append("window horizon must be the canonical 3 s P3 horizon")
    if not math.isclose(float(d.get("sample_period_s", math.nan)), DT_S, rel_tol=0.0, abs_tol=1e-12):
        f.append("sample period must be the canonical 5 ms source period")
    if d.get("complete_window_samples") != SAMPLES:
        f.append(f"complete_window_samples must equal {SAMPLES}")
    for key in _REQUIRED_TRUE_FLAGS:
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in _FORBIDDEN_TRUE_FLAGS:
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    representation = d.get("finite_window_representation")
    if representation not in ("oscillator_shaping_state", "equivalent_hard_finite_window_constraint"):
        f.append("finite-window representation is not an allowed SEA3 theorem representation")
    if d.get("provider_generated_source_family") is not True:
        f.append("provider_generated_source_family is not true")
    if d.get("front_end_entry_witness_id") in (None, ""):
        f.append("missing same-history front-end entry witness")
    if d.get("live_covariance_seed_witness_id") in (None, ""):
        f.append("missing same-history Live covariance seed witness")
    f.extend(_continuity_failures(d.get("transitions")))
    return list(dict.fromkeys(f))


def validate_artifact(d: dict[str, Any]) -> list[str]:
    """Canonical P3 acceptance gate for a finite SEA3 window artifact."""
    f = validate_candidate_structure(d)
    if not PROVIDER_IMPLEMENTATION_CLOSED:
        f.append("validated SEA0 hard finite-window provider is not implemented")
    return list(dict.fromkeys(f))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict[str, Any]:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    physical = PHYSICAL.build(path)
    response = RESPONSE.directional_response_enclosure(REPO)
    prerequisite_failures = {
        "complete": COMPLETE.validate(complete),
        "physical": PHYSICAL.validate(physical),
        "response": RESPONSE.validate(response),
    }
    prerequisite_failures = {k: v for k, v in prerequisite_failures.items() if v}
    if prerequisite_failures:
        raise RuntimeError(f"SEA3 hard-window prerequisites failed: {prerequisite_failures}")

    executable = {
        "machine_readable_R_lambda_closed": MACHINE_READABLE_R_LAMBDA_CLOSED,
        "hard_shaping_state_or_excitation_bound_closed": HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED,
        "joint_translational_rotational_shaping_closed": JOINT_TRANSLATIONAL_ROTATIONAL_SHAPING_CLOSED,
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "SEA3_parameter_domain_compact": True,
        "compact_transition_relation_is_theorem_domain": True,
        "compactness_is_not_an_open_obligation": True,
        "allowed_finite_window_representations": [
            "oscillator_shaping_state",
            "equivalent_hard_finite_window_constraint",
        ],
        "window_horizon_s": HORIZON_S,
        "sample_period_s": DT_S,
        "complete_window_samples": SAMPLES,
        "executable_provider_ingredients": executable,
        "provider_implementation_closed": PROVIDER_IMPLEMENTATION_CLOSED,
        "finite_window_realization_certificate_closed": False,
        "source_reachable_event_family_materialized": False,
        "continuum_RAO_set_theorem_retained": True,
        "continuum_RAO_set_theorem_is_not_finite_window_generator": True,
        "trajectory_replay_used": False,
        "gaussian_good_event_used": False,
        "spectral_moment_only_source_used": False,
        "arbitrary_bounded_input_source_used": False,
        "fixed_lambda_word_used": False,
        "independent_axis_boxes_used": False,
        "independent_SEA_RAO_cartesian_product_used": False,
        "finite_RAO_grid_used": False,
        "independent_tuner_schedule_used": False,
        "retired_P2_graph_used": False,
        "selected_four_S_word_used": False,
        "P3_promoted": False,
        "next_obligation": (
            "machine-execute the already-compact SEA3 theorem domain: encode and validate R_lambda, "
            "a hard x^s shaping-state/excitation or equivalent hard-window constraint, and the joint "
            "translational/rotational response from that same phase-continuous history over all 601 samples"
        ),
    }


def validate_status(d: dict[str, Any]) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("canonical source mismatch")
    for key in (
        "SEA3_parameter_domain_compact",
        "compact_transition_relation_is_theorem_domain",
        "compactness_is_not_an_open_obligation",
        "continuum_RAO_set_theorem_retained",
        "continuum_RAO_set_theorem_is_not_finite_window_generator",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    ingredients = d.get("executable_provider_ingredients", {})
    expected = {
        "machine_readable_R_lambda_closed": MACHINE_READABLE_R_LAMBDA_CLOSED,
        "hard_shaping_state_or_excitation_bound_closed": HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED,
        "joint_translational_rotational_shaping_closed": JOINT_TRANSLATIONAL_ROTATIONAL_SHAPING_CLOSED,
    }
    if ingredients != expected:
        f.append("executable provider ingredient gates drifted")
    if d.get("provider_implementation_closed") is not PROVIDER_IMPLEMENTATION_CLOSED:
        f.append("provider gate mismatch")
    for key in (
        "finite_window_realization_certificate_closed",
        "source_reachable_event_family_materialized",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false before SEA0 finite-window closure")
    for key in _FORBIDDEN_TRUE_FLAGS:
        if d.get(key) is not False:
            f.append(f"status reintroduced forbidden shortcut {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate_status(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "SEA3_compact": d["SEA3_parameter_domain_compact"],
        "compactness_open": not d["compactness_is_not_an_open_obligation"],
        "executable_provider_ingredients": d["executable_provider_ingredients"],
        "provider_closed": d["provider_implementation_closed"],
        "family_materialized": d["source_reachable_event_family_materialized"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
