#!/usr/bin/env python3
"""Non-promoting SEA3 screening estimates for the OU-III stability proof.

This tool intentionally does not construct a stability certificate.  It only
turns shipping source clamps and the WavePeriodEstimator's declared horizons
into deterministic sizing numbers for the SEA0/P2/P3 work.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "OU3_SEA3_INITIAL_ESTIMATES_V1"
M_MAX = 3


def _const_expr(text: str, name: str) -> str:
    match = re.search(
        rf"\b{name}\s*=\s*([^;]+);",
        text,
    )
    if not match:
        raise ValueError(f"could not find {name}")
    return match.group(1).strip()


def _simple_float_expr(expr: str) -> float:
    """Evaluate only a literal or one literal division used by source constants."""
    cleaned = expr.replace("f", "").replace("F", "").strip()
    if "/" in cleaned:
        parts = [part.strip() for part in cleaned.split("/")]
        if len(parts) != 2:
            raise ValueError(f"unsupported expression: {expr}")
        return float(parts[0]) / float(parts[1])
    return float(cleaned)


def _const_float(text: str, name: str) -> float:
    return _simple_float_expr(_const_expr(text, name))


def _wave_period_defaults(text: str) -> tuple[float, float, float]:
    match = re.search(
        r"moment_horizon_periods\s*=\s*([0-9.]+)f.*?"
        r"min_horizon_sec\s*=\s*([0-9.]+)f.*?"
        r"max_horizon_sec\s*=\s*([0-9.]+)f",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise ValueError("could not find WavePeriodEstimator default horizons")
    return tuple(float(value) for value in match.groups())  # type: ignore[return-value]


def build_estimates(repo_root: Path) -> dict[str, Any]:
    filter_header = (
        repo_root / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
    ).read_text(encoding="utf-8")
    limits_header = (
        repo_root / "src/tuner/SeaStateAdaptationLimits.h"
    ).read_text(encoding="utf-8")
    period_header = (
        repo_root / "src/tuner/WavePeriodEstimator.h"
    ).read_text(encoding="utf-8")

    dt = _const_float(filter_header, "FREQ_SMOOTHER_DT")
    f_min = _const_float(filter_header, "MIN_TUNE_FREQ_HZ")
    f_max = _const_float(filter_header, "MAX_TUNE_FREQ_HZ")
    tau_max = _const_float(filter_header, "MAX_TAU_S")
    sigma_max = _const_float(filter_header, "MAX_SIGMA_A")
    rs_min = _const_float(filter_header, "MIN_R_S")
    rs_max = _const_float(filter_header, "MAX_R_S")
    pseudo_max = _const_float(
        filter_header, "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"
    )
    ema_max = _const_float(
        limits_header, "kDynamicEmaHorizonMaxSec"
    )
    moment_periods, moment_floor, moment_ceiling = _wave_period_defaults(
        period_header
    )

    tune_period_min = 1.0 / f_max
    tune_period_max = 1.0 / f_min

    def phase_deg(freq_hz: float, duration_s: float) -> float:
        return round(360.0 * freq_hz * duration_s, 12)

    # Screening proxy only: the exact theorem window will be derived from the
    # three-mode directional sea/RAO enclosure and its finite-window
    # observability certificate.
    screening_window_min = min(
        moment_ceiling,
        max(moment_floor, moment_periods * tune_period_min),
    )
    screening_window_max = min(
        moment_ceiling,
        max(moment_floor, moment_periods * tune_period_max),
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "proof_status": "exploratory_non_promoting",
        "certificate_promoted": False,
        "trajectory_replay_used": False,
        "source_generated": True,
        "m_max": M_MAX,
        "sea_class": "up_to_three_directional_JONSWAP_PM_partitions",
        "p2_relationship": {
            "current_interface": "OU3_P2_CORRELATED_STAGE_TRANSFER_V1",
            "intent": "sound_refinement_not_replacement",
            "required_inclusion": (
                "L_actual_sea subset Lhat_sea3 subset L_current_source"
            ),
        },
        "shipping_clamps": {
            "imu_dt_s": dt,
            "tune_frequency_hz": {"min": f_min, "max": f_max},
            "tau_aw_s": {"max": tau_max},
            "sigma_aw_mps2": {"max": sigma_max},
            "r_s": {"min": rs_min, "max": rs_max},
            "s_pseudo_update_period_s": {"max": pseudo_max},
            "dynamic_ema_horizon_s": {"max": ema_max},
            "normal_live_accelerometer_rejection_allowed": False,
        },
        "derived_screening_estimates": {
            "tune_period_s": {
                "min": round(tune_period_min, 12),
                "max": round(tune_period_max, 12),
            },
            "single_component_phase_advance_per_imu_sample_deg": {
                "min": phase_deg(f_min, dt),
                "max": phase_deg(f_max, dt),
            },
            "single_component_phase_advance_over_max_s_gap_deg": {
                "min": phase_deg(f_min, pseudo_max),
                "max": phase_deg(f_max, pseudo_max),
            },
            "max_s_gap_imu_samples": int(math.ceil(pseudo_max / dt)),
            "wave_period_estimator_screening_proxy": {
                "definition": (
                    "clamp(moment_horizon_periods / f_tune, "
                    "min_horizon_s, max_horizon_s)"
                ),
                "moment_horizon_periods": moment_periods,
                "absolute_horizon_s": {
                    "min": moment_floor,
                    "max": moment_ceiling,
                },
                "screening_window_s": {
                    "min": round(screening_window_min, 12),
                    "max": round(screening_window_max, 12),
                },
                "screening_window_imu_samples": {
                    "min": int(math.ceil(screening_window_min / dt)),
                    "max": int(math.ceil(screening_window_max / dt)),
                },
                "promotion_use": False,
            },
            "oscillator_shaping_state_examples": [
                {
                    "oscillator_pairs_per_mode": q,
                    "modes": M_MAX,
                    "sea_shaping_state_dimension": 2 * q * M_MAX,
                }
                for q in (2, 4, 6)
            ],
        },
        "not_yet_certified": [
            "physical H_s/T_p/gamma/direction/spreading domain",
            "sea-parameter rate bounds",
            "directional vessel/IMU response (RAO) enclosure",
            "hard finite-window sea IQC or equivalent oscillator enclosure",
            "SEA3-to-P2 reachable-language enclosure",
            "P3 finite-window information/covariance margin",
            "P4 nonlinear lifted dissipation margin",
            "P5 finite capture from the 45-degree entrance",
        ],
    }


def _json_text(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed JSON does not match current source clamps",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()

    payload = build_estimates(args.repo_root)
    rendered = _json_text(payload)
    artifact = args.repo_root / "tools/ou3_sea3_initial_estimates.json"

    if args.check:
        if not artifact.exists():
            raise SystemExit(f"missing {artifact}")
        if artifact.read_text(encoding="utf-8") != rendered:
            raise SystemExit(
                "SEA3 initial estimates are stale; rerun "
                "python3 tools/ou3_sea3_initial_estimates.py"
            )
        return 0

    artifact.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
