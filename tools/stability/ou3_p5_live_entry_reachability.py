#!/usr/bin/env python3
"""What the deployed startup path can actually hand to Live, from the code.

The P4 entry set has been carried as an INDEPENDENT product box of declared
startup handoff radii: attitude, gyro bias, 5 m/s velocity, 20 m position,
300 m*s integral displacement, .3 g latent acceleration, .4 m/s^2 bias error,
each free of the others.  The deployed startup path does not produce that set,
and the reason is structural rather than statistical.

Three facts, each checked here against the shipping sources rather than
restated:

1.  ``updateFrontEnd`` runs the whole front end with ``drive_mekf=false``, and
    every MEKF call inside ``updateCore_`` is guarded by that flag.  So
    ``time_update``/``measurement_update_acc_only`` are never invoked before
    ``goLive``: the MEKF translational block, the latent-acceleration state and
    both bias states are never propagated during startup and stay at the
    constructor value, which ``xext.setZero()`` sets to zero.

2.  ``goLive`` calls ``initialize_from_attitude``, which seats the attitude on
    the Mahony proxy quaternion and calls ``zero_AL_cross_cov_once_()``.  That
    zeroes the ENTIRE attitude-block to linear-block cross covariance, so at
    Live entrance ``P_theta,v = P_theta,p = P_theta,S = P_theta,aw = 0`` and in
    particular the first ``S=0`` Joseph event has ``K_theta,S = 0`` exactly.

3.  ``enterLive_`` reseats the latent OU covariance on the committed operating
    point via ``reset_aw_covariance_to_stationary()``.

The consequence is a hard, fully correlated entry relation.  Writing
``e = x_hat - x_true`` at the Live entrance instant ``T``:

    e_v(T)  = -v_true(T),      e_p(T)  = -p_true(T),
    e_S(T)  = -S_true(T),      e_aw(T) = -a_w_true(T),
    e_bg(T) = -b_g_true(T),    e_b(T)  = -b_true(T),

with the attitude error being whatever the Mahony/proxy stage delivers.  Every
non-attitude entry coordinate is the negated physical truth at one instant of
one admitted BRMM history.  They are not independent coordinates at all: they
are one point of one trajectory, and ``e_S`` is the exact running integral of
``e_p`` along that same trajectory.

This does NOT shrink the P4 basin for proof convenience.  It is the
qualification argument the entry set was missing, and it moves the open
question to where it belongs: the reachable entry radii are the BRMM bounded
motion primitives ``V_m``, ``P_m``, ``S_m`` and the BIAS true-bias envelope,
and ``V_m``, ``P_m``, ``S_m`` are declared but NOT instantiated in the physical
BRMM contract.  Instantiating them is the class-E source qualification the
end-to-end theorem needs; inventing values for them here would be exactly the
fitting the contract forbids.

The audit is over the whole fusion source, not one function: every MEKF entry
point that writes the error state is enumerated, and any occurrence outside the
``drive_mekf`` guard other than the three argued above fails validation instead
of being absorbed silently.

Scope.  This is a one-instant reachability lemma about the deployed handoff.
It says nothing about WHEN ``goLive`` fires or how large the proxy attitude
error is: finite-time capture and the attitude radius remain the P5 obligation,
and no part of P5 is assumed by P4 here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2
import ou3_brmm_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
FUSION = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
QUALIFICATION = "OU3_P5_DEPLOYED_LIVE_ENTRY_REACHABILITY_V1"

# Entry coordinates the deployed handoff pins to the negated physical truth.
PINNED_COORDINATES = ("velocity", "position", "integral_displacement",
                      "latent_acceleration", "gyro_bias", "accelerometer_bias")

# BRMM bounded-motion primitives the pinned coordinates need.
REQUIRED_BRMM_PRIMITIVES = ("V_m", "P_m", "S_m")

FUSION_PARITY = {
    "front_end_does_not_drive_mekf":
        "updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);",
    "mekf_time_update_is_guarded": "mekf_->time_update(gyro, dt);",
    "mekf_accelerometer_update_is_guarded":
        "mekf_->measurement_update_acc_only(acc_in, tempC);",
    "go_live_seats_attitude_from_proxy":
        "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);",
    "go_live_enters_live": "enterLive_();",
    "enter_live_reseats_latent_covariance":
        "mekf_->reset_aw_covariance_to_stationary();",
}
# Every MEKF entry point that writes the error state.  A future edit that adds
# one of these to the startup path has to be re-argued, not silently absorbed.
MEKF_STATE_WRITERS = (
    "time_update", "measurement_update_acc_only", "measurement_update_mag_only",
    "applyIntegralZeroPseudoMeas", "measurement_update_position_pseudo",
    "measurement_update_velocity_pseudo", "measurement_update_vert_velocity_pseudo",
    "initialize_from_attitude", "initialize_from_acc", "initialize_from_acc_preserve_yaw",
    "initialize_from_truth", "set_quaternion_boat",
)
# The three state writers that legitimately sit outside the ``drive_mekf`` guard,
# each with its own argument for why it cannot disturb the linear block.
#
#   initialize_from_attitude   IS the Live entrance event, and it zeroes the
#                              attitude-to-linear cross covariance itself.
#   initialize_from_acc        reachable pre-Live through the public passthrough.
#                              It writes ``xext.head<3>()`` only -- the attitude
#                              error -- and also calls zero_AL_cross_cov_once_(),
#                              so it leaves the linear and bias states untouched
#                              and strengthens the cross-covariance claim.
#   measurement_update_mag_only reachable pre-Live through the public updateMag.
#                              The constructor writes only block-diagonal
#                              covariance seeds, so ``P_lin,att = 0``; the only
#                              propagator that could create it is the guarded
#                              ``time_update``; and a Joseph update whose linear
#                              gain rows are ``P_lin,att S^{-1} = 0`` leaves both
#                              the linear states and that cross at zero.
ALLOWED_UNGUARDED_STATE_WRITERS = frozenset(
    {"initialize_from_attitude", "initialize_from_acc", "measurement_update_mag_only"})

MEKF_PARITY = {
    "attitude_init_zeroes_attitude_linear_cross": "zero_AL_cross_cov_once_();",
    "constructor_zeroes_the_whole_covariance": "Pext.setZero();",
    "initial_linear_uncertainty_is_block_diagonal_only":
        "Pext.template block<3,3>(OFF_S, OFF_S) = Matrix3::Identity() * (sigma_S0 * sigma_S0);",
    "cross_zeroing_covers_all_twelve_linear_states":
        "Pext.template block<NA,NL>(0, OFF_V).setZero();",
    "cross_zeroing_is_symmetric":
        "Pext.template block<NL,NA>(OFF_V, 0).setZero();",
    "constructor_zeroes_the_error_state": "xext.setZero();",
    "linear_block_is_twelve_states": "constexpr int NL = 12; // [v,p,S,a_w]",
}


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _guarded_mekf_calls(text: str) -> bool:
    """Every MEKF drive call sits inside an ``if (drive_mekf)`` block.

    The shipping body is small and regular: the guarded region opens with
    ``if (drive_mekf) {`` and the two MEKF drive calls appear only inside such a
    region.  Checking the brace depth relative to the nearest preceding guard is
    enough to reject a future edit that moves one of them out.
    """
    lines = text.splitlines()
    guard_depth = None
    depth = 0
    seen = {"time": False, "acc": False}
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("if (drive_mekf)"):
            guard_depth = depth
        if "mekf_->time_update(gyro, dt);" in line:
            if guard_depth is None or depth <= guard_depth:
                return False
            seen["time"] = True
        if "mekf_->measurement_update_acc_only(acc_in, tempC);" in line:
            if guard_depth is None or depth <= guard_depth:
                return False
            seen["acc"] = True
        depth += line.count("{") - line.count("}")
        if guard_depth is not None and depth <= guard_depth:
            guard_depth = None
    return all(seen.values())


def _unguarded_state_writers(text: str) -> set[str]:
    """MEKF state-writing calls in the whole file that sit outside the guard."""
    found: set[str] = set()
    guard_depth = None
    depth = 0
    for line in text.splitlines():
        if line.strip().startswith("if (drive_mekf)"):
            guard_depth = depth
        for name in MEKF_STATE_WRITERS:
            if "mekf_->" + name + "(" in line:
                if guard_depth is None or depth <= guard_depth:
                    found.add(name)
        depth += line.count("{") - line.count("}")
        if guard_depth is not None and depth <= guard_depth:
            guard_depth = None
    return found


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    fusion_text = FUSION.read_text(encoding="utf-8")
    mekf_text = MEKF.read_text(encoding="utf-8")

    parity = {k: (v in fusion_text) for k, v in FUSION_PARITY.items()}
    parity.update({k: (v in mekf_text) for k, v in MEKF_PARITY.items()})
    parity["mekf_drive_calls_are_inside_the_guard"] = _guarded_mekf_calls(fusion_text)
    unguarded = _unguarded_state_writers(fusion_text)
    parity["no_unexpected_unguarded_state_writer"] = unguarded <= ALLOWED_UNGUARDED_STATE_WRITERS
    parity_pass = all(parity.values())

    b0, b1, b2 = BIAS0.build(), BIAS1.build(), BIAS2.build()
    bad = {k: v for k, v in (("BIAS0", BIAS0.validate(b0)), ("BIAS1", BIAS1.validate(b1)),
                             ("BIAS2", BIAS2.validate(b2))) if v}
    if bad:
        raise RuntimeError("live-entry reachability prerequisites failed: " + repr(bad))
    true_bias_upper = max(float(b["true_bias_norm_upper_mps2"]) for b in (b0, b1, b2))

    source = COMPLETE.build(Path(domain_path))
    src_failures = COMPLETE.validate(source)
    if src_failures:
        raise RuntimeError("complete BRMM source invalid: " + repr(src_failures))
    primitives = source["physical_BRMM_contract"]["unfrozen_physical_constants"]
    missing = [k for k in REQUIRED_BRMM_PRIMITIVES if primitives.get(k) is None]

    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    declared = {
        "velocity": float(handoff["velocity_error_norm_upper_mps"]),
        "position": float(handoff["position_error_norm_upper_m"]),
        "integral_displacement": float(handoff["integral_displacement_error_norm_upper_m_s"]),
        "latent_acceleration": float(handoff["latent_acceleration_error_norm_upper_mps2"]),
        "gyro_bias": float(handoff["gyro_bias_error_norm_upper_rad_s"]),
        "accelerometer_bias": float(handoff["accelerometer_bias_error_norm_upper_mps2"]),
    }
    # What each pinned coordinate is actually equal to at the handoff instant,
    # and which declared source supplies its radius.
    pinned = {
        "velocity": {"equals": "-v_true(T)", "radius_source": "BRMM V_m",
                     "instantiated": primitives.get("V_m") is not None,
                     "value": primitives.get("V_m")},
        "position": {"equals": "-p_true(T)", "radius_source": "BRMM P_m",
                     "instantiated": primitives.get("P_m") is not None,
                     "value": primitives.get("P_m")},
        "integral_displacement": {"equals": "-S_true(T)", "radius_source": "BRMM S_m",
                                  "instantiated": primitives.get("S_m") is not None,
                                  "value": primitives.get("S_m")},
        "latent_acceleration": {
            "equals": "-a_w_true(T)", "radius_source": "declared non-gravitational CoG cap",
            "instantiated": True,
            "value": float(domain["normal_live"]["non_gravitational_cog_acceleration_norm_upper_mps2"])},
        "gyro_bias": {"equals": "-b_g_true(T)", "radius_source": "declared startup gyro bias bound",
                      "instantiated": True,
                      "value": float(domain["startup"]["initial_tangent_gyro_bias_norm_upper_rad_s"])},
        "accelerometer_bias": {"equals": "-b_true(T)", "radius_source": "BIAS0/1/2 true-bias envelope",
                               "instantiated": True, "value": true_bias_upper},
    }
    if set(pinned) != set(PINNED_COORDINATES):
        raise RuntimeError("pinned coordinate table lost a deployed coordinate")

    tightened = {k: bool(v["instantiated"] and float(v["value"]) < declared[k])
                 for k, v in pinned.items()}
    open_radii = sorted(k for k, v in pinned.items() if not v["instantiated"])

    entry_relation = {
        "attitude": "Mahony/proxy attitude error at goLive; the P5 capture obligation",
        "every_other_coordinate": "negated physical truth of ONE admitted BRMM history at ONE instant",
        "coordinates_are_independent": False,
        "integral_is_the_running_integral_of_the_same_position_history": True,
        "cross_covariance_attitude_to_linear_at_entry": 0.0,
        "first_S_zero_event_attitude_gain": 0.0,
    }

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": parity_pass,
        "fusion_sha256": hashlib.sha256(FUSION.read_bytes()).hexdigest(),
        "mekf_sha256": hashlib.sha256(MEKF.read_bytes()).hexdigest(),

        "mekf_state_writers_outside_the_drive_guard": sorted(unguarded),
        "allowed_unguarded_state_writers": sorted(ALLOWED_UNGUARDED_STATE_WRITERS),
        "pre_live_accel_relock_writes_only_the_attitude_head": (
            "initialize_from_acc sets xext.head<3>() to zero and calls "
            "zero_AL_cross_cov_once_(); it never writes the linear or bias states"),
        "pre_live_magnetometer_cannot_move_the_linear_block": (
            "constructor writes only block-diagonal covariance seeds, so P_lin,att = 0; "
            "the only propagator that could create it is the guarded time_update, and a "
            "Joseph update with zero linear gain rows leaves states and cross at zero"),
        "mekf_linear_block_propagated_before_go_live": False,
        "mekf_states_at_go_live": "constructor value, zero, for v/p/S/a_w and both bias blocks",
        "attitude_linear_cross_covariance_zeroed_at_go_live": True,
        "latent_covariance_reseated_at_go_live": True,
        "deployed_live_entry_relation": entry_relation,

        "pinned_entry_coordinates": pinned,
        "declared_independent_handoff_radii": declared,
        "deployed_pin_is_tighter_than_declared": tightened,
        "entry_radii_still_open": open_radii,
        "required_BRMM_primitives": list(REQUIRED_BRMM_PRIMITIVES),
        "uninstantiated_BRMM_primitives": missing,
        "BRMM_motion_primitives_instantiated": not missing,
        "inventing_primitive_values_here": False,

        "reachability_argument_is_from_deployed_code": True,
        "reachability_argument_is_replay_or_fit": False,
        "P5_capture_assumed": False,
        "go_live_timing_established_here": False,
        "proxy_attitude_radius_established_here": False,
        "P4_promoted_here": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "instantiate the BRMM bounded-motion primitives V_m, P_m and S_m from the "
            "declared source class, then prove finite-time Mahony/proxy capture: goLive "
            "fires within a bounded time with a certified proxy attitude radius"),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("source changed")
    extra = set(d.get("mekf_state_writers_outside_the_drive_guard", ())) - ALLOWED_UNGUARDED_STATE_WRITERS
    if extra:
        f.append("unargued MEKF state writer outside the drive guard: " + repr(sorted(extra)))
    if set(d.get("allowed_unguarded_state_writers", ())) != set(ALLOWED_UNGUARDED_STATE_WRITERS):
        f.append("allowed unguarded writer list changed without re-argument")
    if not d.get("pre_live_magnetometer_cannot_move_the_linear_block"):
        f.append("pre-Live magnetometer argument dropped")
    if not d.get("pre_live_accel_relock_writes_only_the_attitude_head"):
        f.append("pre-Live accelerometer re-lock argument dropped")
    if d.get("shipping_source_parity_pass") is not True:
        f.append("shipping source parity lost: " +
                 repr(sorted(k for k, v in d.get("shipping_source_parity", {}).items() if not v)))
    for k in ("attitude_linear_cross_covariance_zeroed_at_go_live",
              "latent_covariance_reseated_at_go_live",
              "reachability_argument_is_from_deployed_code"):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in ("mekf_linear_block_propagated_before_go_live",
              "reachability_argument_is_replay_or_fit", "P5_capture_assumed",
              "go_live_timing_established_here", "proxy_attitude_radius_established_here",
              "inventing_primitive_values_here", "P4_promoted_here", "P5_MAY_START"):
        if d.get(k) is not False:
            f.append(k + " not false")
    rel = d.get("deployed_live_entry_relation", {})
    if rel.get("coordinates_are_independent") is not False:
        f.append("entry relation claims independent coordinates")
    if rel.get("integral_is_the_running_integral_of_the_same_position_history") is not True:
        f.append("integral/position kinematic tie lost")
    if float(rel.get("cross_covariance_attitude_to_linear_at_entry", 1.0)) != 0.0:
        f.append("attitude/linear cross covariance at entry is not zero")
    if float(rel.get("first_S_zero_event_attitude_gain", 1.0)) != 0.0:
        f.append("first S=0 attitude gain is not zero")
    pinned = d.get("pinned_entry_coordinates", {})
    if set(pinned) != set(PINNED_COORDINATES):
        f.append("pinned coordinate table incomplete")
    for name, row in pinned.items():
        if not str(row.get("equals", "")).startswith("-"):
            f.append(name + " is not pinned to a negated truth value")
        if row.get("instantiated") and not (float(row.get("value", -1.0)) > 0.0):
            f.append(name + " radius is not positive")
    # The open radii must be exactly the uninstantiated BRMM primitives.
    open_radii = set(d.get("entry_radii_still_open", ()))
    expected = {n for n, row in pinned.items() if not row.get("instantiated")}
    if open_radii != expected:
        f.append("open entry-radius list does not match the uninstantiated primitives")
    if bool(d.get("BRMM_motion_primitives_instantiated")) != (not d.get("uninstantiated_BRMM_primitives")):
        f.append("primitive instantiation flag inconsistent with its own list")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "parity": d["shipping_source_parity_pass"],
        "linear_propagated_before_go_live": d["mekf_linear_block_propagated_before_go_live"],
        "pinned": {k: v["equals"] for k, v in d["pinned_entry_coordinates"].items()},
        "tighter_than_declared": d["deployed_pin_is_tighter_than_declared"],
        "uninstantiated_primitives": d["uninstantiated_BRMM_primitives"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
