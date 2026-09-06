#!/usr/bin/env python3
"""Non-promoting fixed-history complete-SEA3 source-core feasibility point.

This is a numerical realization of the *existing* continuum Hilbert-ball SEA3
source from :mod:`ou3_sea3_continuum_driver_gram`; it is not another source
model. The purpose is to supply one legal, sustained same-history member to the
ledger's mandatory complete-word feasibility experiment.

One exact admissible SEA3 history is fixed:

* one active JONSWAP partition, H=1.5 m, Tp=6 s, gamma=3.3;
* constant lambda over a 60 s startup prehistory and the 601-sample diagnostic
  window that immediately follows it;
* direction-independent translational RAO
      h_z(f)=G min(1,(fc/f)^2), G=1, fc=0.5 Hz,
  with h_x=h_y=0;
* zero rotational response, an admissible member of the separately declared
  Normal-Live rotational family;
* one common continuum driver for prehistory and diagnostic window.

The hard driver is chosen analytically, not from a replay or a frequency grid:

    a_1(f,theta) = beta sqrt(D(theta)) g(f),
    a_2=a_3=0,
    g(f)=C exp(-(f-f0)^2/(2 sigma_f^2)),
    integral |g(f)|^2 df = 1,
    beta=0.5.

Because the directional density is normalized, integral D(theta)dtheta=1,
||a||_H=beta<=1 exactly. The same a drives all times and channels. With a
direction-independent RAO the theta integral collapses analytically, leaving
only the continuum frequency integral for the vertical CoG acceleration.

The Simpson evaluations below are numerical evaluation of that continuum
integral. Quadrature nodes are never source modes or independently selectable
amplitudes. Two resolutions are compared as a non-promoting numerical
convergence check. Mathematical source membership comes from the analytic
Hilbert norm above.

The first diagnostic sample is evaluated at absolute source time 60 s, exactly
where the startup prehistory ends. There is no phase reset at handoff. This
module deliberately stops before inventing a front-end state or covariance
seed. The next stage runs [0,60 s) of this same member through shipping C++
updateFrontEnd/TunerReady/goLive, serializes the actual Live entry, and executes
the same source from t=60 s through every due S=0 operation at its actual-applied
per-axis R_S.
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
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_FIXED_HISTORY_SOURCE_CORE_V3"
N = 601
DT = 0.005
PREHISTORY_S = 60.0
WINDOW_START_S = PREHISTORY_S
H_M = 1.5
TP_S = 6.0
GAMMA = 3.3
RAO_GAIN = 1.0
RAO_CORNER_HZ = 0.5
RAO_POWER = 2.0
DRIVER_CENTER_HZ = 1.0 / TP_S
DRIVER_SIGMA_HZ = 0.005
DRIVER_BETA = 0.5
SIGMA_LO = 0.07
SIGMA_HI = 0.09
PM_EXPONENT = 1.25


def _shape(f_hz: float) -> float:
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


def _driver_raw(f_hz: float) -> float:
    z = (f_hz - DRIVER_CENTER_HZ) / DRIVER_SIGMA_HZ
    return math.exp(-0.5 * z * z)


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
        value = fn(f) * f
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


def _normalized_spectrum_scale(panels: int) -> float:
    lo, hi = _integration_limits()
    raw = _simpson_log(_shape, lo, hi, panels)
    if not (raw > 0.0 and math.isfinite(raw)):
        raise RuntimeError("JONSWAP normalization integral failed")
    return (H_M * H_M / 16.0) / raw


def _driver_normalization(panels: int) -> float:
    lo, hi = _integration_limits()
    norm2 = _simpson_log(lambda f: _driver_raw(f) ** 2, lo, hi, panels)
    if not (norm2 > 0.0 and math.isfinite(norm2)):
        raise RuntimeError("continuum driver normalization failed")
    return 1.0 / math.sqrt(norm2)


def acceleration_at_time(t_s: float, panels: int = 4096) -> float:
    """Evaluate the exact fixed continuum member at one absolute source time."""
    lo, hi = _integration_limits()
    spectrum_scale = _normalized_spectrum_scale(panels)
    driver_c = _driver_normalization(panels)

    def integrand(f_hz: float) -> float:
        spectral = math.sqrt(max(0.0, spectrum_scale * _shape(f_hz)))
        driver = DRIVER_BETA * driver_c * _driver_raw(f_hz)
        return spectral * _acc_transfer(f_hz) * driver * math.cos(2.0 * math.pi * f_hz * t_s)

    return _simpson_log(integrand, lo, hi, panels)


def _evaluate_window(panels: int) -> list[float]:
    return [acceleration_at_time(WINDOW_START_S + k * DT, panels) for k in range(N)]


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
    RESPONSE.evaluate_rao_envelope_member(RAO_GAIN, RAO_CORNER_HZ, RAO_POWER, response)

    y_coarse = _evaluate_window(1024)
    y_fine = _evaluate_window(2048)
    max_abs_delta = max(abs(a - b) for a, b in zip(y_coarse, y_fine))
    max_abs = max(abs(x) for x in y_fine)
    normal_live_accel_cap = float(
        json.loads(domain_path.read_text(encoding="utf-8"))["normal_live"]
        ["non_gravitational_cog_acceleration_norm_upper_mps2"]
    )

    source_core = [
        {
            "k": k,
            "source_time_s": WINDOW_START_S + k * DT,
            "word_time_s": k * DT,
            "f_cog_body_mps2": [0.0, 0.0, y_fine[k]],
            "omega_body_corrected_rad_s": [0.0, 0.0, 0.0],
        }
        for k in range(N)
    ]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "role": "non-promoting legal fixed-history point for the ledger complete-word feasibility experiment",
        "source_membership": {
            "hard_driver_qualification": driver["qualification"],
            "fixed_history_operator": driver["fixed_history_operator"],
            "driver_field": "a1(f,theta)=beta*sqrt(D(theta))*C*exp(-(f-f0)^2/(2*sigma_f^2)); a2=a3=0",
            "driver_center_hz": DRIVER_CENTER_HZ,
            "driver_sigma_hz": DRIVER_SIGMA_HZ,
            "driver_beta": DRIVER_BETA,
            "driver_norm": DRIVER_BETA,
            "directional_norm_identity": "integral D(theta)dtheta=1",
            "membership_is_analytic_Hilbert_norm_not_quadrature": True,
            "same_driver_field_prehistory_and_window": True,
            "same_driver_field_translation_and_rotation": True,
        },
        "SEA3_fixed_history": {
            "active_partitions": 1,
            "H_r_m": [H_M, 0.0, 0.0],
            "Tp_r_s": [TP_S, TP_S, TP_S],
            "gamma_r": [GAMMA, 1.0, 1.0],
            "lambda_constant_over_prehistory_and_window": True,
            "partition_peak_steepness_admissible": True,
            "total_Hs_m": H_M,
            "directional_density_integrates_to_one": True,
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
        "same_history_phase_contract": {
            "prehistory_start_source_time_s": 0.0,
            "prehistory_end_source_time_s": PREHISTORY_S,
            "word_start_source_time_s": WINDOW_START_S,
            "prehistory_duration_s": PREHISTORY_S,
            "phase_reset_at_live_entry": False,
            "word_uses_immediate_phase_continuation": True,
        },
        "quadrature_diagnostic": {
            "coordinate": "log frequency",
            "source_modes_are_quadrature_nodes": False,
            "coarse_panels": 1024,
            "fine_panels": 2048,
            "max_abs_sample_delta_coarse_to_fine": max_abs_delta,
            "max_abs_source_acceleration_mps2": max_abs,
            "normal_live_acceleration_cap_mps2": normal_live_accel_cap,
            "inside_normal_live_acceleration_cap": max_abs <= normal_live_accel_cap,
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
            "run source times [0,60 s) through shipping C++ updateFrontEnd until/through TunerReady, call the real goLive handoff exactly at source time 60 s, serialize the resulting frontend and H18/A21 covariance state, then execute source time 60 s onward through all 601 samples and every actual-applied R_S event"
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
        "membership_is_analytic_Hilbert_norm_not_quadrature",
        "same_driver_field_prehistory_and_window",
        "same_driver_field_translation_and_rotation",
    ):
        if membership.get(key) is not True:
            f.append(f"source membership lost {key}")
    norm = float(membership.get("driver_norm", math.inf))
    if not (0.0 < norm <= 1.0):
        f.append("fixed-history continuum driver escaped the admitted hard unit ball")
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
    phase = d.get("same_history_phase_contract", {})
    if phase.get("phase_reset_at_live_entry") is not False:
        f.append("fixed source reseeded phase at Live entry")
    if phase.get("word_uses_immediate_phase_continuation") is not True:
        f.append("fixed source word is not the prehistory phase continuation")
    if not math.isclose(float(phase.get("prehistory_duration_s", -1.0)), PREHISTORY_S):
        f.append("fixed source prehistory duration drifted")
    quad = d.get("quadrature_diagnostic", {})
    if quad.get("source_modes_are_quadrature_nodes") is not False:
        f.append("quadrature nodes were promoted to source modes")
    if not float(quad.get("convergence_relative_to_peak", math.inf)) < 2e-4:
        f.append("fixed-history continuum quadrature has not converged sufficiently for feasibility")
    if quad.get("inside_normal_live_acceleration_cap") is not True:
        f.append("selected fixed source point leaves the declared Normal-Live acceleration domain")
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
        "membership": d["source_membership"],
        "history": d["SEA3_fixed_history"],
        "phase": d["same_history_phase_contract"],
        "response": d["fixed_response_member"],
        "quadrature": d["quadrature_diagnostic"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
