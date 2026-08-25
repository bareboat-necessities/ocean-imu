#!/usr/bin/env python3
"""Source-faithful heading branch contract for the OU-III P5 handoff.

P1's normal/timeout cosine bounds are gravity-direction (tilt) bounds.  They are
not full SO(3) attitude bounds.  The shipping startup wrapper has two distinct
heading branches:

* in the full-heading, with-magnetometer quality handoff, ``north_ready`` is
  required and the pending absolute yaw gauge is composed with the proxy tilt;
* the timeout handoff does not require ``north_ready``.  It may therefore enter
  H mode with the proxy's unobservable yaw and the large free-yaw covariance.

Consequently P5 may only build a full-heading Cayley node for a source branch
that actually owns a heading gauge.  The ungauged timeout branch must remain on
the gravity/yaw-quotient theorem route until a magnetic-lock/regauge hybrid event
supplies that gauge.

For a gauged branch, if P1 certifies cos(theta_t)>=c_t>=0 and the declared
heading-gauge error is |psi|<=psi_bar<1, then

  cos(theta_full)
    >= c_t sqrt(1-psi_bar^2) - sqrt(1-c_t^2) psi_bar.

The inequalities cos(psi)>=sqrt(1-psi^2) and sin(psi)<=psi keep the arithmetic
elementary and conservative.  The resulting full-attitude Cayley bound is then
computed from |c|^2=4(1-cos theta)/(1+cos theta).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def sqrt_down(x: float) -> float:
    if x < 0.0:
        raise RuntimeError("negative sqrt")
    y = math.sqrt(x)
    return math.nextafter(y, -math.inf)


def sqrt_up(x: float) -> float:
    if x < 0.0:
        raise RuntimeError("negative sqrt")
    y = math.sqrt(x)
    return math.nextafter(y, math.inf)


def _gauged_full_cosine_lower(tilt_cos_lower: float, psi_upper: float) -> float:
    c = float(tilt_cos_lower)
    p = float(psi_upper)
    if not (0.0 <= c <= 1.0 and 0.0 <= p < 1.0):
        raise RuntimeError("gauged full-attitude composition requires c_t in [0,1], psi in [0,1)")
    cos_psi_lo = sqrt_down(down(1.0 - up(p*p)))
    sin_tilt_hi = sqrt_up(up(max(0.0, 1.0 - down(c*c))))
    return down(down(c * cos_psi_lo) - up(sin_tilt_hi * p))


def _cayley_from_cos_lower(c: float) -> float:
    if not (-1.0 < c <= 1.0):
        raise RuntimeError("full attitude bound reaches Cayley singularity")
    num = up(1.0 - c)
    den = down(1.0 + c)
    return up(2.0 * sqrt_up(up(num / den)))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("heading handoff domain must not be trajectory fitted")
    p1 = P1.build(domain_path)
    failures = [f"P1: {x}" for x in P1.validate(p1)]

    text = WRAPPER.read_text(encoding="utf-8")
    required = {
        "quality_north_gate": "const bool north_ready = !cfg_.with_mag || mag_ref_set_;",
        "quality_uses_north_gate": "tilt_trusted &&\n            north_ready &&\n            impl_.isTunerReady()",
        "timeout_does_not_name_north_gate": "const bool ready_by_timeout =\n            proxy_ready &&\n            (t_ >= timeout_sec) &&\n            mag_gravity_aligned_branch_;",
        "gauge_presence_test": "const bool have_yaw_gauge = std::isfinite(pending_yaw_abs_rad_);",
        "gauged_seed": "? boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)\n                : q_proxy;",
        "free_yaw_sigma_branch": "? cfg_.proxy_handoff_yaw_sigma_rad\n            : cfg_.proxy_handoff_yaw_sigma_free_rad;",
    }
    for label, marker in required.items():
        if marker not in text:
            failures.append(f"missing source heading semantic {label}")

    psi = float(p1["go_live"]["full_heading_internal_gauge_error_upper_rad"])
    cn = float(p1["normal_handoff"]["true_gravity_cosine_lower"])
    ct = float(p1["timeout_handoff"]["combined_true_gravity_cosine_lower"])
    cfn = _gauged_full_cosine_lower(cn, psi)
    cft = _gauged_full_cosine_lower(ct, psi)
    qn = _cayley_from_cos_lower(cfn)
    qt = _cayley_from_cos_lower(cft)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SOURCE_FAITHFUL_HEADING_HANDOFF_BRANCH_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "P1_gravity_cosines_are_tilt_only": True,
        "full_heading_quantitative_scope_requires_magnetometer": True,
        "heading_gauge_error_upper_rad": psi,
        "gauged_quality_handoff": {
            "with_magnetometer_quality_branch_requires_north_ready": True,
            "tilt_cosine_lower": cn,
            "full_attitude_cosine_lower": cfn,
            "full_attitude_cayley_norm_upper": qn,
            "full_heading_P5_node_available": True,
        },
        "gauged_timeout_subbranch": {
            "tilt_cosine_lower": ct,
            "full_attitude_cosine_lower": cft,
            "full_attitude_cayley_norm_upper": qt,
            "full_heading_P5_node_available": True,
        },
        "ungauged_timeout_subbranch": {
            "source_timeout_requires_north_ready": False,
            "source_can_handoff_with_pending_yaw_abs_nan": True,
            "full_heading_cayley_bound_available": False,
            "yaw_covariance_seed_rad": float(p1["go_live"]["ungauged_yaw_covariance_sigma_rad"]),
            "required_route": "GRAVITY_ONLY_YAW_QUOTIENT_UNTIL_MAGNETIC_GAUGE_HYBRID_EVENT",
        },
        "timeout_full_heading_node_covers_complete_timeout_family": False,
        "pass": not failures and qn < 1.0 and qt < 1.0,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("heading handoff is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("heading handoff uses replay")
    if d.get("filter_changed") is not False:
        failures.append("heading handoff changes filter")
    if d.get("P1_gravity_cosines_are_tilt_only") is not True:
        failures.append("P1 gravity cosine was treated as full attitude")
    q = d.get("gauged_quality_handoff", {})
    if q.get("with_magnetometer_quality_branch_requires_north_ready") is not True:
        failures.append("quality full-heading branch lost north gate")
    if q.get("full_heading_P5_node_available") is not True:
        failures.append("gauged quality branch not available")
    u = d.get("ungauged_timeout_subbranch", {})
    if u.get("source_timeout_requires_north_ready") is not False:
        failures.append("timeout incorrectly requires north")
    if u.get("full_heading_cayley_bound_available") is not False:
        failures.append("ungauged timeout was assigned a full-heading Cayley bound")
    if "YAW_QUOTIENT" not in str(u.get("required_route", "")):
        failures.append("ungauged timeout does not route through yaw quotient")
    if d.get("timeout_full_heading_node_covers_complete_timeout_family") is not False:
        failures.append("gauged timeout node incorrectly covers complete timeout family")
    for section in ("gauged_quality_handoff", "gauged_timeout_subbranch"):
        row = d.get(section, {})
        c = row.get("full_attitude_cosine_lower")
        qq = row.get("full_attitude_cayley_norm_upper")
        if not (isinstance(c, (int, float)) and -1.0 < float(c) <= 1.0):
            failures.append(f"{section}: invalid full cosine")
        if not (isinstance(qq, (int, float)) and 0.0 < float(qq) < 1.0):
            failures.append(f"{section}: invalid Cayley node")
    if d.get("pass") is not True:
        failures.append("heading handoff contract did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "quality": out.get("gauged_quality_handoff"),
        "gauged_timeout": out.get("gauged_timeout_subbranch"),
        "ungauged_timeout": out.get("ungauged_timeout_subbranch"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
