#!/usr/bin/env python3
"""V2 adapter for the first P5 rotation-gauged accelerometer stage.

The underlying source family is axis-isotropic in the first S/a_w prefix, but
outward interval arithmetic can widen numerically identical axis expressions by
different ulps.  V1 required the three *interval representations* to compare
bit-identically and therefore stopped before evaluating any correction cell.

V2 keeps the stronger structural requirements that matter to the rotation
gauge: all cross-axis linear covariance entries remain exactly zero and every
theta/a_w entry remains exactly zero.  It merely hulls the three equivalent
axis enclosures for P_SS, P_Saw and P_aw before continuing.  Hulling outward
representations cannot exclude a source value and does not assume equality of
the interval endpoints.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_first_accel_rotation_gauge as V1G

DEFAULT_DOMAIN = V1G.DEFAULT_DOMAIN
SCHEMA = 2


def _scalar_axis_structure_hulled(Pm) -> tuple[Interval, Interval, Interval]:
    pss = hull(*(Pm[12 + ax][12 + ax] for ax in range(3)))
    psa = hull(*(Pm[12 + ax][15 + ax] for ax in range(3)))
    paw = hull(*(Pm[15 + ax][15 + ax] for ax in range(3)))
    for ai in range(3):
        for aj in range(3):
            if ai != aj:
                for a, b in ((12 + ai, 12 + aj), (12 + ai, 15 + aj), (15 + ai, 15 + aj)):
                    z = Pm[a][b]
                    if z.lo != 0.0 or z.hi != 0.0:
                        raise RuntimeError("first-prefix linear covariance gained cross-axis terms")
    for ti in range(3):
        for aj in range(3):
            z = Pm[ti][15 + aj]
            if z.lo != 0.0 or z.hi != 0.0:
                raise RuntimeError("first-prefix theta/a_w covariance is not exactly zero")
    return pss, psa, paw


def _install_backend() -> None:
    V1G._scalar_axis_structure = _scalar_axis_structure_hulled


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    yaw_axis_face_pieces: int = 4,
    force_magnitude_pieces: int = 4,
) -> dict:
    _install_backend()
    out = dict(V1G.build(
        Path(domain_path).resolve(),
        source_pieces=source_pieces,
        yaw_axis_face_pieces=yaw_axis_face_pieces,
        force_magnitude_pieces=force_magnitude_pieces,
    ))
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_FIRST_ACCEL_ROTATION_GAUGED_AXIS_HULLED_SUBDIVISION"
    out["axis_isotropic_source_intervals_hulled_across_equivalent_axes"] = True
    out["bit_identical_axis_interval_endpoints_required"] = False
    out["cross_axis_covariance_still_required_exact_zero"] = True
    out["theta_aw_covariance_still_required_exact_zero"] = True
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V1G.SCHEMA
    failures = V1G.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("axis_isotropic_source_intervals_hulled_across_equivalent_axes") is not True:
        failures.append("axis-equivalent interval hull is not active")
    if d.get("bit_identical_axis_interval_endpoints_required") is not False:
        failures.append("bit-identical interval endpoint gate remains active")
    if d.get("cross_axis_covariance_still_required_exact_zero") is not True:
        failures.append("cross-axis zero structure was weakened")
    if d.get("theta_aw_covariance_still_required_exact_zero") is not True:
        failures.append("theta/a_w zero structure was weakened")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--yaw-axis-face-pieces", type=int, default=4)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain.resolve(),
        source_pieces=args.source_pieces,
        yaw_axis_face_pieces=args.yaw_axis_face_pieces,
        force_magnitude_pieces=args.force_magnitude_pieces,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE"],
        "children": out["evaluated_child_count"],
        "over_limit": out["children_above_validated_correction_limit"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "margin": out["minimum_correction_range_margin_rad"],
        "fixed_inverse": out["fixed_pivot_inverse_count"],
        "fallback_inverse": out["spectral_fallback_inverse_count"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
