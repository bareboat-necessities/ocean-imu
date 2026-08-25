#!/usr/bin/env python3
"""P5 full-H prefix backend with large finite deployed corrections.

V2 removed the broad-tau covariance dependency artifact but still delegated each
attitude correction to the local correction-Cayley helper, whose promoted range
ends at 3 rad.  CI then reached the first accelerometer correction with the
physical error still at only q~=1.036 and stopped solely because the correction
cell exceeded that helper range.

This backend keeps V2's complete 18x18 covariance/Joseph/reset propagation and
replaces only that proof-coordinate operation with
``ou3_p5_deployed_quaternion_cayley_cell``.  The shipping normalized quaternion
is composed directly with the current error quaternion and only the *resulting*
rotation is converted back to Cayley.  Corrections may therefore pass through
pi; the real gate is whether the resulting source cell reaches the antipodal
set.

No estimator code, covariance update, gain, source branch, or theorem domain is
changed.  Numerical nonclosure remains fail-closed and supplies the next source
subdivision witness.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_deployed_quaternion_cayley_cell as QCOMP
import ou3_p5_full_h_prefix_cells as V1
import ou3_p5_full_h_prefix_cells_v2 as V2

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 3


def _install_backend() -> None:
    V2._install_backend()
    V1.SIGNED = QCOMP


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    _install_backend()
    primitive = QCOMP.build(Path(domain_path).resolve())
    pf = QCOMP.validate(primitive)
    out = dict(V1.build(Path(domain_path).resolve()))
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_FULL_18X18_H_PREFIX_DEPLOYED_QUATERNION_PROPAGATION"
    out["active_full_matrix_backend"] = "DEPENDENCY_PRESERVING_OU_KERNEL_BOUNDS_PLUS_DEPLOYED_QUATERNION_COMPOSITION"
    out["integrated_ou_transition_monotone_endpoint_hull_used"] = True
    out["integrated_ou_process_positive_kernel_moment_bounds_used"] = True
    out["goLive_gyro_bias_covariance_includes_full_startup_RW_upper"] = True
    out["correction_cayley_coordinate_formed_before_group_composition"] = False
    out["deployed_quaternion_composed_before_result_cayley"] = True
    out["correction_norm_three_rad_is_promotion_gate"] = False
    out["maximum_validated_deployed_correction_norm_rad"] = primitive["maximum_validated_correction_norm_rad"]
    out["only_resulting_error_antipode_is_group_chart_gate"] = True
    out["deployed_quaternion_primitive_status"] = primitive["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE"]
    out["deployed_quaternion_primitive_failures"] = pf
    if pf:
        out.setdefault("failures", []).extend(f"deployed-quaternion: {x}" for x in pf)
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V1.SCHEMA
    failures = V1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("active_full_matrix_backend") != "DEPENDENCY_PRESERVING_OU_KERNEL_BOUNDS_PLUS_DEPLOYED_QUATERNION_COMPOSITION":
        failures.append("v3 full-matrix backend is not active")
    for k in (
        "integrated_ou_transition_monotone_endpoint_hull_used",
        "integrated_ou_process_positive_kernel_moment_bounds_used",
        "goLive_gyro_bias_covariance_includes_full_startup_RW_upper",
        "deployed_quaternion_composed_before_result_cayley",
        "only_resulting_error_antipode_is_group_chart_gate",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "correction_cayley_coordinate_formed_before_group_composition",
        "correction_norm_three_rad_is_promotion_gate",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if d.get("deployed_quaternion_primitive_status") != "PASS":
        failures.append("deployed quaternion composition prerequisite did not pass")
    failures += list(d.get("deployed_quaternion_primitive_failures", []))
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"],
        "q8_closed": out["complete_q_le_8_prefix_family_closed"],
        "max_q": out["max_reached_cayley_norm_upper"],
        "smaller_chart": out["smaller_source_reachable_chart_upper"],
        "backend": out["active_full_matrix_backend"],
        "max_validated_correction_rad": out["maximum_validated_deployed_correction_norm_rad"],
        "first_failure": out["first_failure"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
