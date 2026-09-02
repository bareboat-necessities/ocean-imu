#!/usr/bin/env python3
"""Probe the whole-word P3 construction on one exact P2 source node.

This is an intermediate numerical certificate, not source-path P3 promotion.
It keeps the existing source-derived covariance upper for the selected node and
replaces only the one-step injected-noise comparison by a word-accumulated
process floor:

* exact integrated-OU Q(Nh) in word scaling for each [v,p,S,a_w] axis;
* exact doubling recursion for each (attitude, gyro-bias) axis; and
* the active accelerometer-bias process floor in A mode.

The source path may change tuner cells inside the word.  Therefore a useful
exact-node result is evidence that the process algebra has enough margin, but
only the clock-phase source-path construction may promote canonical P3.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_scaled_process as SCALED
import ou3_p3_word_process_floor as WORD
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _translation_word_margin(mode, x, sigma, rs, raw, live, vector, process, sched):
    h = float(sched["dt_s"])
    steps_cap = math.floor(BASE.down(float(raw["word_horizon_s_lower"]) / BASE.up(h)))
    if steps_cap < 1:
        raise RuntimeError("word has no certain prediction step")
    max_k = int(math.floor(math.log2(steps_cap)))

    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qba = BASE.pos(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "active BA process",
    )
    cadence_lo = float(raw["cadence_s"][0])
    upper = raw["Sigma_diagonal_upper"]
    physical = [float(upper[6]), float(upper[9]), float(upper[12]), float(upper[15])]

    best = None
    for k in range(max_k, -1, -1):
        steps = 2 ** k
        horizon = BASE.down(steps * h)
        Xlo = BASE.down(steps * x.lo)
        Xhi = BASE.up(steps * x.hi)
        if Xhi > WORD.WORD_EXACT_SERIES_MAX_X:
            continue
        try:
            infos = WORD.translation_information_cover(Interval(Xlo, Xhi))
        except RuntimeError:
            continue

        scales = [
            BASE.down(sigma.lo * horizon),
            BASE.down(sigma.lo * horizon * horizon),
            BASE.down(sigma.lo * horizon * horizon * horizon),
            BASE.down(sigma.lo),
        ]
        if any(v <= 0.0 for v in scales):
            continue
        sigma_root = [
            BASE.up(math.sqrt(BASE.up(physical[i])) / scales[i]) for i in range(4)
        ]

        # Conservative measurement information that can be assimilated over
        # this suffix.  Accelerometer/magnetometer cross information is charged
        # to a_w as a diagonal domination, matching the validated predecessor
        # derivation; nothing is silently discarded.
        s_firings = BASE.up(BASE.up(horizon / BASE.down(cadence_lo)) + 1.0)
        info_S = BASE.up(
            s_firings * BASE.up(scales[2] * scales[2] / BASE.down(rs.lo * rs.lo))
        )
        per_sample_aw = BASE.up(
            (
                fhi * fhi * qtheta
                + sigma.lo * sigma.lo
                + (qba if mode == "A" else 0.0)
            ) / ra
        )
        per_sample_aw = BASE.up(per_sample_aw + BASE.up(mhi * mhi * qtheta / rm))
        info_aw = BASE.up(float(steps) * per_sample_aw)
        scaled_info = [
            0.0,
            0.0,
            BASE.up(info_S * BASE.up(sigma_root[2] * sigma_root[2])),
            BASE.up(info_aw * BASE.up(sigma_root[3] * sigma_root[3])),
        ]

        margin = math.inf
        worst_cell = None
        for cell, information in infos:
            d = WORD.translation_margin_from_information(information, sigma_root, scaled_info)
            if not d > 0.0:
                margin = 0.0
                break
            if d < margin:
                margin = d
                worst_cell = cell.as_list()
        if not margin > 0.0:
            continue
        row = {
            "margin_lower": BASE.down(margin),
            "prediction_steps": steps,
            "horizon_s": horizon,
            "X_horizon_over_tau": [Xlo, Xhi],
            "information_subcells": len(infos),
            "worst_information_subcell": worst_cell,
            "S_firings_upper": s_firings,
        }
        if best is None or row["margin_lower"] > best["margin_lower"]:
            best = row

    if best is None:
        raise RuntimeError("no word suffix certified translation process floor")
    return best


def _attitude_bias_margin(raw, live, vector, process, sched):
    h = float(sched["dt_s"])
    k = WORD.word_step_doublings(float(raw["word_horizon_s_lower"]), h)
    steps = 2 ** k
    ab = process["attitude_gyro_bias"]
    qtheta = BASE.pos(ab["theta_diagonal_lower"], "theta process")
    qbg = BASE.pos(ab["gyro_bias_diagonal_lower"], "bias process")
    cross = float(ab["cross_norm_upper"])
    rho = BASE.down(1.0 - BASE.up(cross / BASE.down(math.sqrt(qtheta * qbg))))
    coupling = BASE.up(h * BASE.up(math.sqrt(BASE.up(qbg / BASE.down(qtheta)))))

    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    per_sample = BASE.up(BASE.up(fhi * fhi / ra + mhi * mhi / rm) * qtheta)

    Omega = WORD.attitude_bias_word_noise(rho, coupling, k, BASE.up(steps * per_sample))
    upper = raw["Sigma_diagonal_upper"]
    Sigma = [
        [WORD.I(BASE.up(upper[0] / qtheta)), WORD.I(0.0)],
        [WORD.I(0.0), WORD.I(BASE.up(upper[3] / qbg))],
    ]
    delta = WORD.generalized_delta(Omega, Sigma, BASE.MIN_USEFUL_DELTA)
    return {
        "margin_lower": delta,
        "prediction_steps": steps,
        "rho_process_lower": rho,
        "coupling_per_step_upper": coupling,
        "attitude_information_per_sample_upper": per_sample,
    }


def _active_ba_margin(mode, raw, sigma, rs, live, vector, process, sched):
    if mode != "A":
        return {"margin_lower": math.inf, "active": False}
    h = float(sched["dt_s"])
    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qba = BASE.pos(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "BA process",
    )
    sS2 = (sigma.lo * h * h * h) ** 2
    betaS = BASE.up(sS2 / BASE.rs_variance_lower(rs, sched))
    betaAcc = BASE.up((fhi * fhi * qtheta + sigma.lo * sigma.lo + qba) / ra)
    betaMag = BASE.up(mhi * mhi * qtheta / rm)
    beta = BASE.up(BASE.up(betaS + betaAcc) + betaMag)
    floor = BASE.down(1.0 / BASE.up(1.0 + beta))
    sigma_upper = BASE.up(
        raw["Sigma_diagonal_upper"][18] / raw["comparison_scale_diagonal_squared"][18]
    )
    return {
        "margin_lower": BASE.down(floor / sigma_upper),
        "active": True,
        "measurement_information_upper": beta,
    }


def _evaluate_mode(mode, xparts, sigma, rs, live, vector, process, sched, alpha6):
    rows = []
    for x, rho in xparts:
        raw = BASE.mode_cell(mode, x, rho, sigma, rs, live, vector, process, sched, alpha6)
        trans = _translation_word_margin(mode, x, sigma, rs, raw, live, vector, process, sched)
        att = _attitude_bias_margin(raw, live, vector, process, sched)
        ba = _active_ba_margin(mode, raw, sigma, rs, live, vector, process, sched)
        delta = BASE.down(min(trans["margin_lower"], att["margin_lower"], ba["margin_lower"]))
        rows.append({
            "x_h_over_tau": x.as_list(),
            "translation": trans,
            "attitude_gyro_bias": att,
            "active_accelerometer_bias": ba,
            "relative_word_process_margin_lower": delta,
        })
    if not rows:
        raise RuntimeError("node generated no validated x subcells")
    worst = min(rows, key=lambda r: r["relative_word_process_margin_lower"])
    return {
        "dimension": 18 if mode == "H" else 21,
        "validated_x_subcell_count": len(rows),
        "relative_word_process_margin_lower": worst["relative_word_process_margin_lower"],
        "useful_margin_gate": BASE.MIN_USEFUL_DELTA,
        "useful_margin_pass": worst["relative_word_process_margin_lower"] >= BASE.MIN_USEFUL_DELTA,
        "worst_subcell": worst,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, source_node_index: int = 0):
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P3 node word probe must not be trajectory fitted")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    node = NODES.node(source_node_index, nodes)

    live = domain["normal_live"]
    vector = BASE.VECTOR.build()
    process = BASE.PROCESS.build()
    words = BASE.WORDS.build(path)
    prereq = []
    prereq += [f"vector: {x}" for x in BASE.VECTOR.validate(vector)]
    prereq += [f"process: {x}" for x in BASE.PROCESS.validate(process)]
    prereq += [f"word: {x}" for x in BASE.WORDS.validate(words)]
    if prereq:
        raise RuntimeError(f"P3 node word prerequisites failed: {prereq}")

    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, node["tau_s"])
    x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    rs = Interval(*map(float, node["R_S_filter_std"]))
    xparts = SCALED.split_x_cell(x)
    alpha6 = BASE.vector_alpha6(live, vector)

    modes = {
        mode: _evaluate_mode(mode, xparts, sigma, rs, live, vector, process, sched, alpha6)
        for mode in ("H", "A")
    }
    failures = []
    for mode in ("H", "A"):
        d = modes[mode]["relative_word_process_margin_lower"]
        if not (math.isfinite(d) and d > 0.0):
            failures.append(f"{mode}: word process margin is not strict")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_EXACT_P2_NODE_WHOLE_WORD_PROCESS_PROBE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_node_index": int(source_node_index),
        "source_node": node,
        "x_h_over_tau": x.as_list(),
        "modes": modes,
        "P3_NODE_WORD_NUMERICAL_PROBE_PASS": not failures,
        "P3_SOURCE_PATH_CERTIFICATE_ESTABLISHED_HERE": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "replace the static-node covariance upper and same-parameter Q(Nh) lower by a clock-phase source-path comparison; "
            "only that pathwise minimum may promote canonical P3"
        ),
        "failures": failures,
    }


def validate(d):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "P3_NODE_WORD_NUMERICAL_PROBE_PASS"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "P3_SOURCE_PATH_CERTIFICATE_ESTABLISHED_HERE", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("relative_word_process_margin_lower")
        if not isinstance(delta, (int, float)) or not math.isfinite(float(delta)) or float(delta) <= 0.0:
            f.append(f"{mode}: invalid word margin")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "node": d["source_node_index"],
        "x": d["x_h_over_tau"],
        "H_delta": d["modes"]["H"]["relative_word_process_margin_lower"],
        "H_useful": d["modes"]["H"]["useful_margin_pass"],
        "A_delta": d["modes"]["A"]["relative_word_process_margin_lower"],
        "A_useful": d["modes"]["A"]["useful_margin_pass"],
        "failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
