#!/usr/bin/env python3
"""Extract the OU-III continuous-source proof domain from implementation guards.

This producer intentionally does not infer bounds from the eight reference
trajectories.  It parses the shipping implementation's safety/clamp constants
and records every discrete branch the validated backend must cover.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HEADER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"

REQUIRED = (
    "MIN_TUNE_FREQ_HZ", "MAX_TUNE_FREQ_HZ", "MIN_TAU_S", "MAX_TAU_S",
    "MAX_SIGMA_A", "MIN_R_S", "MAX_R_S",
    "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT", "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
    "MAG_DELAY_SEC", "ONLINE_TUNE_WARMUP_SEC",
)


def parse_const(text: str, name: str) -> float:
    # Keep the extractor intentionally narrow: a proof-domain change should
    # fail loudly if the implementation stops spelling these as scalar constants.
    pat = re.compile(
        rf"constexpr\s+float\s+{re.escape(name)}\s*=\s*([0-9.+\-eE]+)f?\s*;"
    )
    m = pat.search(text)
    if not m:
        raise RuntimeError(f"cannot extract implementation constant {name}")
    return float(m.group(1))


def build(header: Path) -> dict:
    text = header.read_text()
    c = {name: parse_const(text, name) for name in REQUIRED}
    return {
        "schema": 1,
        "claim": "OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "source_complete_parameter_domain": True,
        "validated_arithmetic": False,
        "outward_rounded": False,
        "implementation_header": str(header.relative_to(REPO)),
        "continuous_parameters": {
            "wave_tune_frequency_hz": [c["MIN_TUNE_FREQ_HZ"], c["MAX_TUNE_FREQ_HZ"]],
            "tau_aw_s": [c["MIN_TAU_S"], c["MAX_TAU_S"]],
            # The implementation has only an upper safety clamp on sigma.  Zero
            # is therefore the source-complete lower endpoint for validation.
            "sigma_aw_mps2": [0.0, c["MAX_SIGMA_A"]],
            "R_S_base": [c["MIN_R_S"], c["MAX_R_S"]],
            "pseudo_update_period_s": [
                c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"],
                c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"],
            ],
        },
        "timing_constants_s": {
            "mag_delay": c["MAG_DELAY_SEC"],
            "online_tune_warmup": c["ONLINE_TUNE_WARMUP_SEC"],
        },
        "discrete_source_branches": {
            "mode": ["H", "A"],
            "accelerometer_gate": ["accepted", "rejected"],
            "magnetometer_gate": ["not_due", "accepted", "rejected"],
            "S_zero_pseudo": ["not_due", "due"],
            "magnetic_gauge": ["unlocked", "locked", "refined"],
            "aw_covariance_sync": ["not_due", "due_psd_increment"],
        },
        "hybrid_obligations": [
            "startup_handoff", "held_to_active", "magnetic_regauge",
            "tilt_reset", "cooldown", "periodic_aw_covariance_sync",
        ],
        "promotion_rule": (
            "these implementation bounds define the source domain only; theorem promotion "
            "still requires outward-rounded interval/Taylor-model propagation over the full domain"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = build(args.header.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
