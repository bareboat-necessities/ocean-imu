#!/usr/bin/env python3
"""Canonical BRMM moving-Riccati tube with stable Q and endpoint covariance.

The numerical integrated-OU backend is the #489 factored implementation, kept
verbatim in ``ou3_brmm_riccati_tube_factored_backend``.  This canonical facade
adds only the endpoint-reference theorem for the global translation covariance
ceiling.  It does not create a second proof route: the same BASE build, source
invariant, current-source interval cover, H/A aggregation and 1e-18 gate are
used.

The endpoint covariance envelope and the optional scalar moving-Riccati
injection margin are deliberately validated separately.  P4 magnitude-only
consumers may use ``validate_covariance_ceiling``; doing so does not assert that
the independent 1e-18 moving-tube contraction gate passed.
"""
from __future__ import annotations

import math

import ou3_brmm_endpoint_covariance as ENDPOINT
import ou3_brmm_riccati_tube_factored_backend as BACKEND


# The backend saved the original BASE.build before installing its stable
# integrated-OU primitives.  Rebinding this one BASE helper therefore changes
# only the covariance reference-time calculation used by that same build.
ENDPOINT.install(BACKEND.BASE)

SCHEMA = BACKEND.SCHEMA
QUALIFICATION = BACKEND.QUALIFICATION
USEFUL_GATE = BACKEND.USEFUL_GATE
DEFAULT_DOMAIN = BACKEND.DEFAULT_DOMAIN

# Preserve the public numerical helpers used by focused tests/debugging.
step_scaled_q = BACKEND.step_scaled_q
step_scaled_q_over_x = BACKEND.step_scaled_q_over_x
split_x_cell = BACKEND.split_x_cell


def build(domain_path=DEFAULT_DOMAIN):
    # Install at call time as well as import time.  Several proof adapters
    # temporarily replace BASE helpers while sharing the same Python module;
    # the canonical producer must reassert the endpoint-reference theorem
    # immediately before evaluating BACKEND's saved original BASE.build.
    ENDPOINT.install(BACKEND.BASE)
    d = BACKEND.build(domain_path)
    timing = d.get("translation_covariance_ceiling", {})
    d["endpoint_referenced_translation_covariance"] = True
    d["post_reconstruction_forward_propagation_used"] = False
    d["covariance_ceiling_argument"] = (
        "finite-memory recurrent vector/S estimator referenced at the word endpoint; "
        "global BRMM adaptive invariant for every nuisance/process upper"
    )
    profile = dict(d.get("numerical_profile", {}))
    profile["endpoint_referenced_translation_covariance"] = True
    profile["post_reconstruction_forward_propagation_removed"] = True
    d["numerical_profile"] = profile
    if timing.get("translation_reference") != "word_endpoint":
        raise RuntimeError("canonical covariance ceiling lost endpoint reference")
    return d


# BASE.main resolves its module-global build at call time.  Point it at this
# facade so CLI and imported use execute exactly the same canonical producer.
BACKEND.BASE.build = build


def validate_covariance_ceiling(payload):
    """Validate only the rigorous endpoint-referenced Pbar envelope.

    This scoped contract intentionally ignores ``moving_Riccati_tube_pass`` and
    the scalar relative injection margin.  Those fields remain governed by the
    full ``validate`` function and are never promoted by this helper.
    """
    failures = []
    if payload.get("schema") != SCHEMA or payload.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "BRMM_dynamic_source_consumed",
        "current_source_interval_cover_only", "time_varying_source_allowed_inside_covariance_memory_window",
        "endpoint_referenced_translation_covariance",
    ):
        if payload.get(key) is not True:
            failures.append(key + " is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "P2_800_state_partition_consumed", "post_reconstruction_forward_propagation_used",
        "P3_PROMOTED", "P4_PROMOTED",
    ):
        if payload.get(key) is not False:
            failures.append(key + " is not false")
    timing = payload.get("translation_covariance_ceiling", {})
    if timing.get("translation_reference") != "word_endpoint":
        failures.append("covariance envelope is not referenced at word endpoint")
    if timing.get("endpoint_referenced_observability") is not True:
        failures.append("endpoint S-observability flag missing")
    if timing.get("endpoint_p_sign_similarity_applied") is not True:
        failures.append("endpoint p-sign similarity flag missing")
    if timing.get("forward_propagation_after_endpoint_reconstruction") is not False:
        failures.append("post-reconstruction forward propagation reappeared")
    cover = payload.get("interval_cover", {})
    for key in ("h_over_tau_leaf_count", "sigma_cells", "R_S_cells", "combined_current_cells"):
        if int(cover.get(key, 0)) <= 0:
            failures.append("invalid covariance-envelope interval cover " + key)
    if int(cover.get("max_split_depth", -1)) != BACKEND.MAX_X_SPLIT_DEPTH:
        failures.append("x split depth metadata mismatch")
    ceiling = payload.get("full_state_covariance_ceiling", {})
    expected = {"H": 18, "A": 21}
    for mode, n in expected.items():
        row = payload.get("modes", {}).get(mode, {})
        p = row.get("Pbar_diagonal_variance_upper")
        c = ceiling.get(mode)
        if not isinstance(p, list) or len(p) != n:
            failures.append(mode + " Pbar dimension invalid")
            continue
        if not isinstance(c, list) or len(c) != n:
            failures.append(mode + " full-state ceiling dimension invalid")
            continue
        for i, (x, y) in enumerate(zip(p, c)):
            xf = float(x); yf = float(y)
            if not (math.isfinite(xf) and xf > 0.0 and math.isfinite(yf) and yf > 0.0):
                failures.append(f"{mode} covariance ceiling entry {i} invalid")
            if xf != yf:
                failures.append(f"{mode} mode Pbar detached from full-state ceiling at {i}")
    return list(dict.fromkeys(failures))


def validate(payload):
    failures = list(BACKEND.validate(payload))
    timing = payload.get("translation_covariance_ceiling", {})
    if payload.get("endpoint_referenced_translation_covariance") is not True:
        failures.append("translation covariance is not endpoint referenced")
    if payload.get("post_reconstruction_forward_propagation_used") is not False:
        failures.append("endpoint covariance was propagated forward a second time")
    if timing.get("translation_reference") != "word_endpoint":
        failures.append("covariance-memory translation reference is not word endpoint")
    if timing.get("endpoint_referenced_observability") is not True:
        failures.append("endpoint S-observability flag missing")
    if timing.get("endpoint_p_sign_similarity_applied") is not True:
        failures.append("endpoint p-sign similarity flag missing")
    if timing.get("forward_propagation_after_endpoint_reconstruction") is not False:
        failures.append("post-reconstruction forward propagation reappeared")
    return list(dict.fromkeys(failures))


def main():
    return BACKEND.main()


if __name__ == "__main__":
    raise SystemExit(main())
