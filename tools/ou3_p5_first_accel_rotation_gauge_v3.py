#!/usr/bin/env python3
"""V3 source-sparsity adapter for the first P5 rotation-gauged accel stage.

The first rotation-gauged implementation correctly identified an exact source
symmetry, but its generic interval 18x18 multiplication widens algebraic zero
products by one IEEE-754 subnormal on every operation.  Those tiny intervals
are not reachable cross-axis covariance; they are natural-extension arithmetic
noise around the exact structural value zero.

V3 does not introduce a numerical zero threshold.  Before propagating any
child it verifies directly on the active V3 source matrices that:

* the initial linear covariance is block diagonal by physical axis;
* the linear transition never mixes physical axes;
* the process covariance has no cross-axis linear entries;
* the complete base/linear initial, transition and process cross blocks are
  exactly zero.

Those source identities imply exact zero cross-axis and theta/a_w covariance
through the first prediction.  V3 may therefore intersect the generic interval
result with those exact structural identities.  Same-axis P_SS, P_Saw and P_aw
intervals are still hulled across the equivalent axes; no source value is
removed and no replay data or fitted tolerance is used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_first_accel_rotation_gauge as G1
import ou3_p5_full_h_prefix_cells as FULL1
import ou3_p5_full_h_prefix_cells_v3 as FULL3

DEFAULT_DOMAIN = G1.DEFAULT_DOMAIN
SCHEMA = 3

_SOURCE_SPARSITY_CERTIFIED = False
_SOURCE_CELLS_CHECKED = 0
_MAX_ARITHMETIC_ZERO_DUST = 0.0


def _exact_zero(x: Interval) -> bool:
    return x.lo == 0.0 and x.hi == 0.0


def _require_zero(x: Interval, label: str) -> None:
    if not _exact_zero(x):
        raise RuntimeError(f"source sparsity lost exact zero at {label}: {x.as_list()}")


def _assert_first_prefix_source_sparsity(domain_path: Path, source_pieces: int) -> int:
    """Prove exact axis/base-linear decoupling before generic interval products."""
    FULL3._install_backend()
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    checked = 0
    for src, phase in G1._source_phase_children(source_pieces):
        P0 = FULL1._initial_covariance(src, domain_path)
        F, Q, _ = FULL1._transition_and_Q(src, domain)

        axis_ids = [[6 + a, 9 + a, 12 + a, 15 + a] for a in range(3)]
        for a in range(3):
            for b in range(3):
                if a == b:
                    continue
                for i in axis_ids[a]:
                    for j in axis_ids[b]:
                        _require_zero(P0[i][j], f"P0 axis {a}->{b} {i},{j}")
                        _require_zero(Q[i][j], f"Q axis {a}->{b} {i},{j}")
                        _require_zero(F[i][j], f"F axis {a}->{b} {i},{j}")

        # goLive decouples the complete attitude/gyro-bias block from the
        # complete linear block.  The first prediction preserves that because
        # F and Q are block diagonal across the same partition.
        for i in range(6):
            for j in range(6, 18):
                _require_zero(P0[i][j], f"P0 base-linear {i},{j}")
                _require_zero(P0[j][i], f"P0 linear-base {j},{i}")
                _require_zero(Q[i][j], f"Q base-linear {i},{j}")
                _require_zero(Q[j][i], f"Q linear-base {j},{i}")
                _require_zero(F[i][j], f"F base-linear {i},{j}")
                _require_zero(F[j][i], f"F linear-base {j},{i}")
        checked += 1
    if checked == 0:
        raise RuntimeError("source sparsity proof evaluated no first-prefix cells")
    return checked


def _scalar_axis_structure_from_exact_sparsity(Pm) -> tuple[Interval, Interval, Interval]:
    global _MAX_ARITHMETIC_ZERO_DUST
    if not _SOURCE_SPARSITY_CERTIFIED:
        raise RuntimeError("first-prefix structural zeros used before source sparsity certificate")

    pss = hull(*(Pm[12 + ax][12 + ax] for ax in range(3)))
    psa = hull(*(Pm[12 + ax][15 + ax] for ax in range(3)))
    paw = hull(*(Pm[15 + ax][15 + ax] for ax in range(3)))

    # The generic natural extension may surround the exact structural zero by
    # subnormal dust.  Require that it still contains zero, record its size for
    # audit, and intersect it with the exact source identity rather than using a
    # floating numerical threshold.
    for ai in range(3):
        for aj in range(3):
            if ai == aj:
                continue
            for a, b in ((12 + ai, 12 + aj), (12 + ai, 15 + aj), (15 + ai, 15 + aj)):
                z = Pm[a][b]
                if not (z.lo <= 0.0 <= z.hi):
                    raise RuntimeError("generic cross-axis enclosure no longer contains structural zero")
                _MAX_ARITHMETIC_ZERO_DUST = max(_MAX_ARITHMETIC_ZERO_DUST, z.abs_upper())
    for ti in range(3):
        for aj in range(3):
            z = Pm[ti][15 + aj]
            if not (z.lo <= 0.0 <= z.hi):
                raise RuntimeError("generic theta/a_w enclosure no longer contains structural zero")
            _MAX_ARITHMETIC_ZERO_DUST = max(_MAX_ARITHMETIC_ZERO_DUST, z.abs_upper())
    return pss, psa, paw


def _install_backend(domain_path: Path, source_pieces: int) -> None:
    global _SOURCE_SPARSITY_CERTIFIED, _SOURCE_CELLS_CHECKED, _MAX_ARITHMETIC_ZERO_DUST
    _SOURCE_SPARSITY_CERTIFIED = False
    _MAX_ARITHMETIC_ZERO_DUST = 0.0
    _SOURCE_CELLS_CHECKED = _assert_first_prefix_source_sparsity(domain_path, source_pieces)
    _SOURCE_SPARSITY_CERTIFIED = True
    G1._scalar_axis_structure = _scalar_axis_structure_from_exact_sparsity


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    yaw_axis_face_pieces: int = 4,
    force_magnitude_pieces: int = 4,
) -> dict:
    domain_path = Path(domain_path).resolve()
    _install_backend(domain_path, source_pieces)
    out = dict(G1.build(
        domain_path,
        source_pieces=source_pieces,
        yaw_axis_face_pieces=yaw_axis_face_pieces,
        force_magnitude_pieces=force_magnitude_pieces,
    ))
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_FIRST_ACCEL_ROTATION_GAUGED_SOURCE_SPARSITY_SUBDIVISION"
    out["first_prefix_source_sparsity_certified_before_interval_product"] = True
    out["source_sparsity_cells_checked"] = _SOURCE_CELLS_CHECKED
    out["structural_zero_canonicalization_uses_numeric_threshold"] = False
    out["cross_axis_interval_dust_treated_as_physical_covariance"] = False
    out["axis_equivalent_same_axis_intervals_hulled"] = True
    out["arithmetic_zero_dust_abs_upper_max"] = _MAX_ARITHMETIC_ZERO_DUST
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = G1.SCHEMA
    failures = G1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("first_prefix_source_sparsity_certified_before_interval_product") is not True:
        failures.append("source sparsity was not certified before canonicalization")
    if int(d.get("source_sparsity_cells_checked", 0)) <= 0:
        failures.append("source sparsity certificate checked no cells")
    if d.get("structural_zero_canonicalization_uses_numeric_threshold") is not False:
        failures.append("structural zero canonicalization uses a numeric threshold")
    if d.get("cross_axis_interval_dust_treated_as_physical_covariance") is not False:
        failures.append("arithmetic zero dust was treated as reachable covariance")
    if d.get("axis_equivalent_same_axis_intervals_hulled") is not True:
        failures.append("same-axis equivalent intervals were not hulled")
    dust = d.get("arithmetic_zero_dust_abs_upper_max")
    if not isinstance(dust, (int, float)) or not math.isfinite(float(dust)) or float(dust) < 0.0:
        failures.append("arithmetic zero dust diagnostic is invalid")
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
        "source_sparsity_cells": out["source_sparsity_cells_checked"],
        "zero_dust_max": out["arithmetic_zero_dust_abs_upper_max"],
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
