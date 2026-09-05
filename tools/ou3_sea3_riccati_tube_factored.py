#!/usr/bin/env python3
"""Canonical SEA3 moving-Riccati tube with stable Q and endpoint covariance.

The numerical integrated-OU backend is the #489 factored implementation, kept
verbatim in ``ou3_sea3_riccati_tube_factored_backend``.  This canonical facade
adds only the endpoint-reference theorem for the global translation covariance
ceiling.  It does not create a second proof route: the same BASE build, source
invariant, current-source interval cover, H/A aggregation and 1e-18 gate are
used.
"""
from __future__ import annotations

import ou3_sea3_endpoint_covariance as ENDPOINT
import ou3_sea3_riccati_tube_factored_backend as BACKEND


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
    d = BACKEND.build(domain_path)
    timing = d.get("covariance_memory", {})
    d["endpoint_referenced_translation_covariance"] = True
    d["post_reconstruction_forward_propagation_used"] = False
    d["covariance_ceiling_argument"] = (
        "finite-memory recurrent vector/S estimator referenced at the word endpoint; "
        "global SEA3 adaptive invariant for every nuisance/process upper"
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


def validate(payload):
    failures = list(BACKEND.validate(payload))
    timing = payload.get("covariance_memory", {})
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
