#!/usr/bin/env python3
"""Non-promoting fixed-history complete-SEA3 source-core feasibility point.

This is the first numerical realization used by the complete-word feasibility
experiment.  It is *not* a second source model.  The source remains the
continuum Hilbert-ball representation from ``ou3_sea3_continuum_driver_gram``.

We choose one exact admissible member of the existing compact SEA3 family:

* one active JONSWAP partition, H=1.5 m, Tp=6 s, gamma=3.3;
* constant lambda over the diagnostic history;
* a direction-independent translational RAO
      h_z(f)=G min(1,(fc/f)^2), G=1, fc=0.5 Hz,
  with h_x=h_y=0;
* zero rotational response, which is admissible under the separately declared
  Normal-Live rotational envelope;
* the same continuum hard driver for the whole window.

For this fixed history let K be the exact continuum map from normalized driver
to the 601 sampled vertical CoG accelerations and Q=K K*.  The chosen driver is

    a = K* e_0 / ||K* e_0||,

so ||a||=1 and the resulting exact source sequence is

    y = Q e_0 / sqrt(Q_00).

Therefore source membership follows from the hard-driver operator itself, not
from the numerical quadrature below.  The quadrature merely evaluates the
continuum correlation integrals for this already-defined member.  It is never
interpreted as a finite harmonic source: no quadrature node is a source mode,
phase coordinate, or independently selectable amplitude.

This module stops at the six-coordinate physical source core.  It does not
invent a front-end entry or covariance seed and does not promote P4/P5.  The
next step is to run a sufficiently long same-history prehistory through the
actual C++ startup/front-end path, then feed the resulting Live entry and these
same continuum-driver samples into the connected 601-sample executor.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable

import ou3_sea3_continuum_driver_gram as DRIVER
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_directional_response_family as RESPONSE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_FIXED_HISTORY_SOURCE_CORE_V1"
N = 601
DT = 0.005
H_M = 1.5
TP_S = 6.0
GAMMA = 3.3
RAO_GAIN = 1.0
RAO_CORNER_HZ = 0.5
RAO_POWER = 2.0
SIGMA_LO = 0.07
SIGMA_HI = 0.09
PM_EXPONENT = 1.25


def _shape(f_hz: float) -> float:
    """Unnormalized JONSWAP shape in frequency coordinates."""
    if not (f_hz > 0.0):
        return 0.0
    fp = 1.0 / TP_S
    x = f_hz / fp
    sigma = SIGMA_LO if x <= 1.0 else SIGMA_HI
    peak = math.exp(-((x - 1.0) ** 2) / (2.0 * sigma * sigma))
    return f_hz**-5 * math.exp(-PM_EXPONENT * x**-4) * GAMMA**peak


def _rao(f_hz: float) -> float:
    if f_hz <= RAO_CORNER_HZ:
        return RAO_GAIN
    return RAO_GAIN * (RAO_CORNER_HZ / f_hz) ** RAO_POWER


def _acc_transfer(f_hz: float) -> float:
    omega = 2.0 * math.pi * f_hz
    return -(omega * omega) * _rao(f_hz)


def _simpson_log(fn: Callable[[float], float], lo: float, hi: float, panels: int) -> float:
    if not (0.0 < lo < hi):
        raise ValueError("positive ordered integration interval required")
    if panels <= 0 or panels % 2:
        raise ValueError("Simpson panel count must be positive and even")
    a = math.log(lo)
    b = math.log(hi)
    h = (b - a) / panels
    total = 0.0
    for i in range(panels + 1):
        t = a + i * h
        f = math.exp(t)
        value = fn(f) * f  # df=f dt
        if i == 0 or i == panels:
            w = 1.0
        elif i % 2:
            w = 4.0
        else:
            w = 2.0
        total += w * value
    return total * h / 3.0


def _integration_limits() -> tuple[float, float]:
    fp = 1.0 / TP_S
    return fp / 64.0, fp * 256.0


def _normalized_spectrum_scale(panels: int) -> tuple[float, float]:
    lo, hi = _integration_limits()
    raw = _simpson_log(_shape, lo, hi, panels)
    if not (raw > 0.0 and math.isfinite(raw)):
        raise RuntimeError("JONSWAP normalization integral failed")
    target_m0 = H_M * H_M / 16.0
    return target_m0 / raw, raw


def _gram_lag(lag: int, panels: int, scale: float) -> float:
    lo, hi = _integration_limits()
    tau = float(lag) * DT
    return _simpson_log(
        lambda f: scale * _shape(f) * (_acc_transfer(f) ** 2) * math.cos(2.0 * math.pi * f * tau),
        lo,
        hi,
        panels,
    )


def _evaluate(panels: int) -> tuple[list[float], dict]:
    scale, raw = _normalized_spectrum_scale(panels)
    q = [_gram_lag(k, panels, scale) for k in range(N)]
    if not (q[0] > 0.0 and math.isfinite(q[0])):
        raise RuntimeError("fixed-history Gram lost positive variance")
    root = math.sqrt(q[0])
    y = [x / root for x in q]
    return y, {
        "panels": panels,
        "raw_shape_integral": raw,
        "spectrum_scale": scale,
        "Q00_acceleration_variance": q[0],
        "driver_norm": 1.0,
        "construction": "a=K*e0/sqrt(e0^T Q e0); y=Qe0/sqrt(Q00)",
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    physical = PHYSICAL.build(domain_path)
    response = RESPONSE.directional_response_enclosure(REPO)
    driver = DRIVER.build()
    bad = {
        "physical": PHYSICAL.validate(physical),
        "response": RESPONSE.validate(response),
        "driver": DRIVER.validate(driver),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"fixed SEA3 source prerequisites failed: {bad}")

    g = float(physical["gravity_mps2"])
    if not PHYSICAL.partition_admissible(H_M, TP_S, g):
        raise RuntimeError("selected fixed SEA3 partition violates peak-steepness contract")
    RESPONSE.evaluate_rao_envelope_member(
        RAO_GAIN, RAO_CORNER_HZ, RAO_POWER, response
    )

    # Two continuum quadrature resolutions.  Agreement is a numerical
    # feasibility diagnostic only; mathematical membership came from y=Ka above.
    y_coarse, coarse = _evaluate(4096)
    y_fine, fine = _evaluate(8192)
    max_abs_delta = max(abs(a - b) for a, b in zip(y_coarse, y_fine))
    max_abs = max(abs(x) for x in y_fine)

    # Identity attitude, zero rotation and a vertical-only legal CoG response.
    source_core = [
        {
            "k": k,
            "t_s": k * DT,
            "f_cog_body_mps2": [0.0, 0.0, y_fine[k]],
            "omega_body_corrected_rad_s": [0.0, 0.0, 0.0],
        }
        for k in range(N)
    ]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "role": "non-promoting legal fixed-history source-core point for the ledger feasibility experiment",
        "source_membership": {
            "hard_driver_qualification": driver["qualification"],
            "fixed_history_operator": driver["fixed_history_operator"],
            "fixed_history_gram": driver["fixed_history_gram"],
            "driver_choice": "a=K*e0/||K*e0||",
            "driver_norm": 1.0,
            "membership_is_by_operator_construction_not_quadrature": True,
            "same_driver_field_entire_window": True,
            "same_driver_field_translation_and_rotation": True,
        },
        "SEA3_fixed_history": {
            "active_partitions": 1,
            "H_r_m": [H_M, 0.0, 0.0],
            "Tp_r_s": [TP_S, TP_S, TP_S],
            "gamma_r": [GAMMA, 1.0, 1.0],
            "lambda_constant_over_window": True,
            "partition_peak_steepness_admissible": True,
            "total_Hs_m": H_M,
            "directional_density_integrates_to_one_and_cancels_for_direction_independent_RAO": True,
        },
        "fixed_response_member": {
            "translation": "h=[0,0,G min(1,(fc/f)^2)]",
            "G": RAO_GAIN,
            "fc_hz": RAO_CORNER_HZ,
            "power": RAO_POWER,
            "rotation": "zero",
            "inside_declared_continuum_RAO_family": True,
            "finite_RAO_grid_used": False,
        },
        "quadrature_diagnostic": {
            "coordinate": "log frequency",
            "source_modes_are_quadrature_nodes": False,
            "coarse": coarse,
            "fine": fine,
            "max_abs_sample_delta_coarse_to_fine": max_abs_delta,
            "max_abs_source_acceleration_mps2": max_abs,
            "convergence_relative_to_peak": max_abs_delta / max(max_abs, 1e-30),
        },
        "sample_count": N,
        "dt_s": DT,
        "source_core": source_core,
        "front_end_entry_derived_from_same_history": False,
        "live_covariance_seed_derived_from_same_history": False,
        "complete_executor_artifact_materialized": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "trajectory_replay_used": False,
        "finite_harmonic_source_used": False,
        "independent_sample_boxes_used": False,
        "next_obligation": (
            "extend this exact fixed continuum-driver member backward through a same-history prehistory, run the actual C++ updateFrontEnd/TunerReady/goLive path, serialize the resulting frontend state and H18/A21 covariance seed, then execute all 601 shipping samples with every actual-applied R_S update"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("fixed history detached from canonical complete SEA3")
    membership = d.get("source_membership", {})
    for key in (
        "membership_is_by_operator_construction_not_quadrature",
        "same_driver_field_entire_window",
        "same_driver_field_translation_and_rotation",
    ):
        if membership.get(key) is not True:
            f.append(f"source membership lost {key}")
    if float(membership.get("driver_norm", 0.0)) != 1.0:
        f.append("fixed-history driver is not on the admitted hard unit ball")
    sea = d.get("SEA3_fixed_history", {})
    if sea.get("active_partitions") != 1:
        f.append("fixed feasibility history changed partition count")
    if sea.get("partition_peak_steepness_admissible") is not True:
        f.append("fixed feasibility history is not physically admissible")
    response = d.get("fixed_response_member", {})
    if response.get("inside_declared_continuum_RAO_family") is not True:
        f.append("fixed response is outside SEA3 RAO family")
    if response.get("finite_RAO_grid_used") is not False:
        f.append("finite RAO grid re-entered fixed-history source")
    quad = d.get("quadrature_diagnostic", {})
    if quad.get("source_modes_are_quadrature_nodes") is not False:
        f.append("quadrature nodes were promoted to source modes")
    if not float(quad.get("convergence_relative_to_peak", math.inf)) < 5e-5:
        f.append("fixed-history continuum quadrature has not converged sufficiently for feasibility")
    if d.get("sample_count") != N or not math.isclose(float(d.get("dt_s", 0.0)), DT):
        f.append("fixed history does not cover canonical 601 samples at 5 ms")
    if not isinstance(d.get("source_core"), list) or len(d["source_core"]) != N:
        f.append("source core payload length mismatch")
    for key in (
        "front_end_entry_derived_from_same_history",
        "live_covariance_seed_derived_from_same_history",
        "complete_executor_artifact_materialized",
        "P4_promoted",
        "P5_promoted",
        "trajectory_replay_used",
        "finite_harmonic_source_used",
        "independent_sample_boxes_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false at source-core stage")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "source": d["canonical_source"],
        "history": d["SEA3_fixed_history"],
        "response": d["fixed_response_member"],
        "quadrature": d["quadrature_diagnostic"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
