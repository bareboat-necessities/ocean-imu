#!/usr/bin/env python3
"""Complete-SEA3 finite process-memory lemma for canonical OU-III P3.

This module does not create a source word and does not use the one-IMU-step
full-state minimum as a contraction ratio.  It extracts one fact that is true
inside every admitted complete SEA3 Normal-Live word: the *active* OU schedule
is piecewise constant because candidate tuning is staged only after
``ADAPT_EVERY_SECS`` and committed at the next IMU boundary.

A conservative integer number of prediction steps therefore exists between two
possible commits.  Over that interval the shipping integrated-OU recursion
injects exactly the finite-horizon covariance Q(H) in each [v,p,S,a_w] axis.
The stable validated integrated-OU backend is evaluated over the whole current
SEA3 tau invariant; sigma uses only its source-uniform lower endpoint.  The
result is then expressed in the same dimensionless translation coordinates used
by the actual-applied-R_S four-S information lemma,

    z = [S, g p, g^2 v, g^3 a_w].

The attitude/gyro-bias and active accelerometer-bias process blocks are treated
over the same physical interval from shipping process densities.  These are
process-memory primitives only.  Measurements, resets and the prior-free
completion are separate same-word operations; this module cannot promote P3 and
cannot replace the SpectralMSE R_S information regularizer.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_four_s_translation_information as FOUR_S
import ou3_sea3_riccati_tube_factored as STABLE_Q
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_FINITE_PROCESS_MEMORY"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    four = FOUR_S.build(path)
    process = PROCESS.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "four_S": FOUR_S.validate(four),
        "process": PROCESS.validate(process),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"complete-SEA3 process-memory prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("process-memory lemma detached from complete SEA3")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    adapt = float(SOURCE.parse_const(wrapper, "ADAPT_EVERY_SECS"))
    dt = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    if not (adapt > 0.0 and dt > 0.0 and adapt > 4.0 * dt):
        raise RuntimeError("shipping adaptation cadence no longer leaves finite process memory")

    parity = {
        "candidate_staged_after_adaptation_cadence": (
            "time_ - last_adapt_time_sec_ > adapt_every_secs_" in wrapper
            and "online_tune_apply_pending_ = true" in wrapper
        ),
        "candidate_committed_at_next_imu_boundary": (
            "void apply_pending_online_tune_()" in wrapper
            and "apply_ou_tune_(false);" in wrapper
            and "online_tune_apply_pending_ = false;" in wrapper
        ),
        "adaptation_cadence_shipping_constant": adapt > 0.0,
    }
    parity_failures = [k for k, v in parity.items() if not v]
    if parity_failures:
        raise RuntimeError(f"shipping staged-commit parity failed: {parity_failures}")

    # Do not depend on exact decimal equality at the > cadence boundary.  One
    # full sample of padding below floor(adapt/dt) leaves an interval that is
    # certainly inside a constant-active-schedule segment.
    nominal_steps = int(math.floor(down(adapt / up(dt))))
    constant_steps = nominal_steps - 1
    if constant_steps < 4:
        raise RuntimeError("constant-active process chunk lost enough samples")
    H = down(constant_steps * dt)

    inv = dynamic["dynamic_invariant"]
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_lo, sigma_hi = map(float, inv["sigma_aw_filter_mps2"])
    if not (0.0 < tau_lo <= tau_hi and 0.0 < sigma_lo <= sigma_hi):
        raise RuntimeError("invalid complete-SEA3 OU invariant")

    x = Interval.outward_bounds(down(H / tau_hi), up(H / tau_lo))
    leaves = STABLE_Q.split_x_cell(x)
    if not leaves:
        raise RuntimeError("stable integrated-OU process cover is empty")
    rho = down(min(float(r) for _, r in leaves))
    if not (math.isfinite(rho) and rho > 0.0):
        raise RuntimeError("finite-memory scaled integrated-OU lower is not strict")

    g = float(four["uniform_S_gap_s_upper"])
    # STABLE_Q uses D=diag(sigma*H, sigma*H^2, sigma*H^3, sigma)
    # in [v,p,S,a].  Transform to [S,g p,g^2 v,g^3 a].
    scale_sq = [
        (sigma_lo * H ** 3) ** 2,
        (g * sigma_lo * H ** 2) ** 2,
        (g ** 2 * sigma_lo * H) ** 2,
        (g ** 3 * sigma_lo) ** 2,
    ]
    q_translation_z = down(rho * min(scale_sq))

    qg = float(process["attitude_gyro_bias"]["q_gyro_lower"])
    qb = float(process["source_constants"]["gyro_bias_rw_variance_density"])
    q_att_bg = down(min(qg * H, qb * H) - up(qb * H * H / 2.0))

    qba_density = float(process["source_constants"]["accel_bias_process_variance_density"])
    tau_ba = float(process["source_constants"]["accel_bias_tau_s"])
    xb = Interval.outward_bounds(down(2.0 * H / tau_ba), up(2.0 * H / tau_ba))
    em1 = VT.expm1_interval(-xb)
    qd_scale = Interval.outward_bounds(-0.5 * tau_ba, -0.5 * tau_ba) * em1
    q_ba = down(qba_density * qd_scale.lo)

    q_H = down(min(q_translation_z, q_att_bg))
    q_A = down(min(q_H, q_ba))
    passed = all(math.isfinite(v) and v > 0.0 for v in (
        rho, q_translation_z, q_att_bg, q_ba, q_H, q_A
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "component_of_complete_SEA3_full_word": True,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_box_used_as_source": False,
        "same_SEA3_dynamic_invariant_consumed": True,
        "shipping_staged_commit_parity": parity,
        "shipping_staged_commit_parity_failures": parity_failures,
        "adapt_every_s": adapt,
        "imu_dt_s": dt,
        "nominal_samples_per_adaptation_cadence": nominal_steps,
        "guaranteed_constant_active_prediction_steps": constant_steps,
        "guaranteed_constant_active_process_horizon_s": H,
        "tau_applied_s": [tau_lo, tau_hi],
        "sigma_aw_filter_mps2": [sigma_lo, sigma_hi],
        "integrated_OU_x_interval": x.as_list(),
        "integrated_OU_validated_leaf_count": len(leaves),
        "integrated_OU_scaled_lambda_min_lower": rho,
        "translation_coordinates": ["S", "g*p", "g^2*v", "g^3*a_w"],
        "translation_coordinate_scale_squared_lower": scale_sq,
        "translation_Q_memory_lambda_min_lower": q_translation_z,
        "attitude_gyro_bias_Q_memory_lambda_min_lower": q_att_bg,
        "active_accel_bias_Q_memory_lambda_min_lower": q_ba,
        "modes": {
            "H": {"dimension": 18, "pre_measurement_process_memory_lambda_min_lower": q_H},
            "A": {"dimension": 21, "pre_measurement_process_memory_lambda_min_lower": q_A},
        },
        "one_step_full_state_Q_used_as_contraction_ratio": False,
        "actual_R_S_information_replaced_by_process_strictness": False,
        "actual_applied_SpectralMSE_R_S_must_be_composed_separately": True,
        "process_memory_pass": passed,
        "P3_promoted": False,
        "next_obligation": (
            "carry this same-word finite process memory through the intervening shipping Joseph/reset "
            "events and combine it with the actual-applied-R_S H18 information matrix in the "
            "prior-free full-matrix completion; then add A-mode finite-bias detectability"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "same_SEA3_dynamic_invariant_consumed",
        "actual_applied_SpectralMSE_R_S_must_be_composed_separately",
        "process_memory_pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_family_replaced", "trajectory_replay_used",
        "independent_tau_sigma_RS_box_used_as_source",
        "one_step_full_state_Q_used_as_contraction_ratio",
        "actual_R_S_information_replaced_by_process_strictness", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("shipping_staged_commit_parity_failures") != []:
        f.append("shipping staged-commit parity failed")
    if int(d.get("guaranteed_constant_active_prediction_steps", 0)) < 4:
        f.append("constant-active process memory too short")
    for mode in ("H", "A"):
        q = d.get("modes", {}).get(mode, {}).get("pre_measurement_process_memory_lambda_min_lower")
        if not isinstance(q, (int, float)) or not (math.isfinite(float(q)) and float(q) > 0.0):
            f.append(f"{mode} process-memory lower is not strict")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "constant_steps": d["guaranteed_constant_active_prediction_steps"],
        "memory_horizon_s": d["guaranteed_constant_active_process_horizon_s"],
        "x_interval": d["integrated_OU_x_interval"],
        "scaled_rho": d["integrated_OU_scaled_lambda_min_lower"],
        "translation_Q_memory": d["translation_Q_memory_lambda_min_lower"],
        "attitude_bg_Q_memory": d["attitude_gyro_bias_Q_memory_lambda_min_lower"],
        "A_bias_Q_memory": d["active_accel_bias_Q_memory_lambda_min_lower"],
        "H_Q_memory": d["modes"]["H"]["pre_measurement_process_memory_lambda_min_lower"],
        "A_Q_memory": d["modes"]["A"]["pre_measurement_process_memory_lambda_min_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
