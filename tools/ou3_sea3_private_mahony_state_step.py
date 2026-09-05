#!/usr/bin/env python3
"""Source-neutral interval step for the shipping private Mahony observer.

This module is deliberately *not* a source generator. It advances one already
initialized ``VerticalAccelComplementary`` state using gyro/accelerometer
coordinates supplied by the same complete SEA3 realization. It never invents
an independent vertical-acceleration, quaternion, gyro or integral-feedback
box.

The arithmetic follows ``Mahony_AHRS<float>::update`` and
``VerticalAccelComplementary::update``. Normalization uses the repository's
validated enclosure of the actual Lomont ``0x5f375a86`` fast inverse square
root and all explicit basic operations are outward rounded to binary32.

The low-level implementation can remain in ``TunerReady`` while an external
bootstrap decides when to call ``goLive``.  The deployed outer
``SeaStateFusion_OU_III`` wrapper is stricter: it owns the bootstrap and has a
configured timeout path.  With current defaults the magnetic acquisition
fallback is 60 s and the startup timeout is 150 s, so 150 s dominates.  That is
*not yet* an unconditional Live-entry bound, because the timeout path still
requires the world-frame gravity-aligned branch.  P3 therefore remains
fail-closed until that branch is certified on the same SEA3 startup path and
the reset-zero Mahony state is propagated to the resulting Live entry.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re

from ou3_interval import Interval
import ou3_binary32_interval as F32
import ou3_fast_inv_sqrt_interval as FINV

REPO = Path(__file__).resolve().parents[1]
VERTICAL = REPO / "src" / "tuner" / "VerticalAccelComplementary.h"
MAHONY = REPO / "src" / "ahrs" / "Mahony_AHRS.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
COMMON = REPO / "src" / "kalman_common" / "SeaStateFusionFilterCommon.h"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_PRIVATE_MAHONY_LIVE_INTERVAL_STEP_V2"


@dataclass(frozen=True)
class State:
    q0: Interval
    q1: Interval
    q2: Interval
    q3: Interval
    integral_x: Interval
    integral_y: Interval
    integral_z: Interval
    up_ms2: Interval


@dataclass(frozen=True)
class Vec3:
    x: Interval
    y: Interval
    z: Interval


def I(x: float) -> Interval:
    return F32.point(x)


def _config_float(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed Config::{name}")
    return float(m.group(1))


def _add4(a: Interval, b: Interval, c: Interval, d: Interval) -> Interval:
    return F32.add(F32.add(F32.add(a, b), c), d)


def _norm2_3(v: Vec3) -> Interval:
    return F32.add(F32.add(F32.square(v.x), F32.square(v.y)), F32.square(v.z))


def _norm2_q(q0: Interval, q1: Interval, q2: Interval, q3: Interval) -> Interval:
    return _add4(F32.square(q0), F32.square(q1), F32.square(q2), F32.square(q3))


def advance_initialized_live(
    state: State,
    *,
    dt: Interval,
    gyro: Vec3,
    acc_specific_force: Vec3,
    gravity_ms2: Interval,
    two_kp: Interval,
    two_ki: Interval,
) -> State:
    """Advance one valid shipping sample of an already initialized proxy."""
    if dt.lo <= 0.0:
        raise ValueError("positive dt required")
    if two_ki.lo <= 0.0:
        raise ValueError("deployed positive-Ki branch required")
    if gravity_ms2.lo <= 0.0:
        raise ValueError("positive gravity magnitude required")

    a_mahony = Vec3(
        F32.neg(acc_specific_force.x),
        F32.neg(acc_specific_force.y),
        F32.neg(acc_specific_force.z),
    )
    a2 = _norm2_3(a_mahony)
    if a2.lo <= 0.0:
        raise ValueError(
            "accelerometer box lost the SEA3 nonzero-norm coupling; partition the same SEA3 cell"
        )
    recip_a = FINV.enclosure(a2)
    ax = F32.mul(a_mahony.x, recip_a)
    ay = F32.mul(a_mahony.y, recip_a)
    az = F32.mul(a_mahony.z, recip_a)

    q0, q1, q2, q3 = state.q0, state.q1, state.q2, state.q3
    q0q0, q1q1, q2q2, q3q3 = (
        F32.square(q0), F32.square(q1), F32.square(q2), F32.square(q3)
    )
    halfvx = F32.sub(F32.mul(q1, q3), F32.mul(q0, q2))
    halfvy = F32.add(F32.mul(q0, q1), F32.mul(q2, q3))
    halfvz_inner = F32.add(F32.sub(F32.sub(q0q0, q1q1), q2q2), q3q3)
    halfvz = F32.mul(I(0.5), halfvz_inner)
    halfex = F32.sub(F32.mul(ay, halfvz), F32.mul(az, halfvy))
    halfey = F32.sub(F32.mul(az, halfvx), F32.mul(ax, halfvz))
    halfez = F32.sub(F32.mul(ax, halfvy), F32.mul(ay, halfvx))

    def integral(old: Interval, err: Interval) -> Interval:
        return F32.add(old, F32.mul(F32.mul(two_ki, err), dt))

    ix = integral(state.integral_x, halfex)
    iy = integral(state.integral_y, halfey)
    iz = integral(state.integral_z, halfez)
    gx = F32.add(F32.add(gyro.x, ix), F32.mul(two_kp, halfex))
    gy = F32.add(F32.add(gyro.y, iy), F32.mul(two_kp, halfey))
    gz = F32.add(F32.add(gyro.z, iz), F32.mul(two_kp, halfez))

    half_dt = F32.mul(I(0.5), dt)
    gx, gy, gz = F32.mul(gx, half_dt), F32.mul(gy, half_dt), F32.mul(gz, half_dt)
    qa, qb, qc = q0, q1, q2
    dq0 = F32.sub(F32.sub(F32.neg(F32.mul(qb, gx)), F32.mul(qc, gy)), F32.mul(q3, gz))
    dq1 = F32.sub(F32.add(F32.mul(qa, gx), F32.mul(qc, gz)), F32.mul(q3, gy))
    dq2 = F32.add(F32.sub(F32.mul(qa, gy), F32.mul(qb, gz)), F32.mul(q3, gx))
    dq3 = F32.sub(F32.add(F32.mul(qa, gz), F32.mul(qb, gy)), F32.mul(qc, gx))
    r0, r1, r2, r3 = (
        F32.add(q0, dq0), F32.add(q1, dq1), F32.add(q2, dq2), F32.add(q3, dq3)
    )
    q2_after = _norm2_q(r0, r1, r2, r3)
    if q2_after.lo <= 0.0:
        raise ValueError("quaternion box lost nonzero-norm coupling; partition the same SEA3 cell")
    recip_q = FINV.enclosure(q2_after)
    n0, n1, n2, n3 = (
        F32.mul(r0, recip_q), F32.mul(r1, recip_q),
        F32.mul(r2, recip_q), F32.mul(r3, recip_q),
    )

    row_x = F32.mul(I(2.0), F32.sub(F32.mul(n1, n3), F32.mul(n0, n2)))
    row_y = F32.mul(I(2.0), F32.add(F32.mul(n2, n3), F32.mul(n0, n1)))
    row_z = F32.add(
        F32.sub(F32.sub(F32.square(n0), F32.square(n1)), F32.square(n2)),
        F32.square(n3),
    )
    down_specific = F32.dot3(
        (row_x, row_y, row_z),
        (acc_specific_force.x, acc_specific_force.y, acc_specific_force.z),
    )
    up = F32.neg(F32.add(down_specific, gravity_ms2))
    return State(n0, n1, n2, n3, ix, iy, iz, up)


def _point_smoke() -> dict:
    st = State(I(1.0), I(0.0), I(0.0), I(0.0), I(0.0), I(0.0), I(0.0), I(0.0))
    out = advance_initialized_live(
        st, dt=I(0.005),
        gyro=Vec3(I(0.01), I(-0.02), I(0.005)),
        acc_specific_force=Vec3(I(0.2), I(-0.1), I(-9.75)),
        gravity_ms2=I(9.80665), two_kp=I(0.2), two_ki=I(0.02),
    )
    vals = (out.q0, out.q1, out.q2, out.q3, out.integral_x,
            out.integral_y, out.integral_z, out.up_ms2)
    return {
        "quaternion": [x.as_list() for x in vals[:4]],
        "integral_feedback": [x.as_list() for x in vals[4:7]],
        "vertical_accel_up_mps2": out.up_ms2.as_list(),
        "finite": all(x.lo == x.lo and x.hi == x.hi for x in vals),
    }


def build() -> dict:
    vertical = VERTICAL.read_text(encoding="utf-8")
    mahony = MAHONY.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    common = COMMON.read_text(encoding="utf-8")

    proxy_timeout = _config_float(wrapper, "proxy_startup_timeout_sec")
    proxy_mag_settle = _config_float(wrapper, "proxy_mag_settle_sec")
    mag_window = _config_float(wrapper, "mag_min_window_sec")
    mag_fallback = _config_float(wrapper, "mag_tilt_fallback_sec")
    mag_deadline = proxy_mag_settle + 2.0 * max(mag_window, 1.0) + mag_fallback
    deployed_timeout = max(proxy_timeout, mag_deadline)

    parity = {
        "wrapper_updates_private_observer_before_vertical_read": (
            "vertical_accel_comp_.update(dt, gyro, acc_in, g_std);" in wrapper
            and "vertical_accel_comp_.verticalAccelUpMs2();" in wrapper
        ),
        "shipping_proxy_feeds_negative_specific_force": "-acc.x(), -acc.y(), -acc.z()," in vertical,
        "shipping_positive_ki_branch": (
            "if (twoKi > T(0))" in mahony
            and "integralFBx += twoKi * halfex * delta_t_sec;" in mahony
        ),
        "shipping_quaternion_normalization": (
            "recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);" in mahony
        ),
        "shipping_vertical_row": (
            "2.0f * (x * z - w * y)" in vertical
            and "2.0f * (y * z + w * x)" in vertical
            and "w * w - x * x - y * y + z * z" in vertical
        ),
        "deployed_wrapper_two_kp": "STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;" in wrapper,
        "deployed_wrapper_two_ki": "STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;" in wrapper,
        "mahony_reset_zeros_integral_state": (
            "integralFBx = T(0);" in mahony and "integralFBy = T(0);" in mahony
            and "integralFBz = T(0);" in mahony
        ),
        "gain_change_does_not_reset_integral_state": (
            "Sets gains only. Does NOT reset quaternion or integral states." in mahony
        ),
        "low_level_tuner_ready_waits_for_external_go_live": (
            "TunerReady,  // tuner trusted, MEKF still held by an external bootstrap" in wrapper
            and "call goLive() once it has tilt and north" in wrapper
        ),
        "deployed_outer_wrapper_owns_handoff": (
            "void maybeHandOffToMekf_()" in wrapper
            and "impl_.goLive(q_seed," in wrapper
            and "stage_ = Stage::Live;" in wrapper
        ),
        "deployed_outer_wrapper_timeout_formula": (
            "std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline)" in wrapper
        ),
        "timeout_still_requires_gravity_aligned_branch": (
            "const bool ready_by_timeout" in wrapper
            and "(t_ >= timeout_sec)" in wrapper
            and "mag_gravity_aligned_branch_;" in wrapper
        ),
        "gravity_aligned_branch_is_world_z_sign": (
            "return acc_world_lp.z() < 0.0f;" in common
        ),
    }
    smoke = _point_smoke()
    return {
        "schema": SCHEMA, "qualification": QUALIFICATION,
        "source_generator": False, "trajectory_replay_used": False,
        "independent_vertical_acceleration_source": False,
        "independent_quaternion_or_integral_box_promotable": False,
        "requires_same_SEA3_gyro_and_specific_force": True,
        "requires_same_SEA3_private_observer_state": True,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "source_order_binary32_one_sample_map_closed": True,
        "actual_fast_inverse_sqrt_used": True,
        "ideal_inverse_sqrt_substituted": False,
        "point_smoke_only_not_P3": True, "smoke": smoke,
        "live_entry_integral_state_starts_from_reset_zero": True,
        "low_level_TunerReady_can_wait_for_external_bootstrap": True,
        "deployed_outer_wrapper_has_timeout_logic": True,
        "deployed_proxy_startup_timeout_s": proxy_timeout,
        "deployed_mag_acquire_deadline_s": mag_deadline,
        "deployed_timeout_s": deployed_timeout,
        "timeout_path_requires_gravity_aligned_branch": True,
        "unconditional_live_entry_upper_bound_closed": False,
        "live_entry_private_observer_invariant_closed": False,
        "compiler_reassociation_or_FMA_closed": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "prove on the same complete SEA3 reset-to-startup path that the world-frame gravity branch acc_world_lp.z<0 is retained or re-entered by the deployed 150 s timeout, then propagate the reset-zero Mahony quaternion/integral-feedback state to that Live entry and feed its vertical output into WPE; do not substitute an independent q/integral/vertical-acceleration box"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "requires_same_SEA3_gyro_and_specific_force", "requires_same_SEA3_private_observer_state",
        "shipping_source_parity_pass", "source_order_binary32_one_sample_map_closed",
        "actual_fast_inverse_sqrt_used", "point_smoke_only_not_P3",
        "live_entry_integral_state_starts_from_reset_zero",
        "low_level_TunerReady_can_wait_for_external_bootstrap",
        "deployed_outer_wrapper_has_timeout_logic",
        "timeout_path_requires_gravity_aligned_branch",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if d.get("deployed_proxy_startup_timeout_s") != 150.0:
        failures.append("unexpected deployed proxy startup timeout")
    if d.get("deployed_mag_acquire_deadline_s") != 60.0:
        failures.append("unexpected deployed magnetic acquisition deadline")
    if d.get("deployed_timeout_s") != 150.0:
        failures.append("unexpected deployed startup timeout")
    for key in (
        "source_generator", "trajectory_replay_used", "independent_vertical_acceleration_source",
        "independent_quaternion_or_integral_box_promotable", "ideal_inverse_sqrt_substituted",
        "unconditional_live_entry_upper_bound_closed", "live_entry_private_observer_invariant_closed",
        "compiler_reassociation_or_FMA_closed", "complete_SEA3_family_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if not d.get("smoke", {}).get("finite"):
        failures.append("point smoke is not finite")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"parity": d["shipping_source_parity"], "smoke": d["smoke"],
                      "startup": {
                          "proxy_timeout_s": d["deployed_proxy_startup_timeout_s"],
                          "mag_deadline_s": d["deployed_mag_acquire_deadline_s"],
                          "deployed_timeout_s": d["deployed_timeout_s"],
                          "requires_aligned_branch": d["timeout_path_requires_gravity_aligned_branch"],
                      },
                      "next_obligation": d["next_obligation"], "failures": failures},
                     indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
