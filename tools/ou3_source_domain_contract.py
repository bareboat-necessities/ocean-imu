#!/usr/bin/env python3
"""Extract the OU-III continuous-source proof domain from implementation guards.

This producer intentionally does not infer bounds from the eight reference
trajectories. It parses the shipping implementation's safety/clamp constants
and records every discrete branch and hybrid transition the validated backend
must cover.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
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

HYBRID_OBLIGATIONS = (
    "startup_handoff",
    "held_to_active",
    "magnetic_lock",
    "magnetic_regauge_refinement",
    "tilt_reset",
    "tilt_relock",
    "cooldown_reentry",
    "periodic_aw_covariance_sync",
)

CONST_RE = re.compile(
    r"constexpr\s+float\s+([A-Za-z_]\w*)\s*=\s*([^;]+);", re.MULTILINE
)


def _strip_float_suffixes(expr: str) -> str:
    return re.sub(r"(?<=\d)[fF]\b", "", expr)


def _eval_constexpr(node: ast.AST, text: str, stack: tuple[str, ...]) -> float:
    if isinstance(node, ast.Expression):
        return _eval_constexpr(node.body, text, stack)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.Name):
        return parse_const(text, node.id, stack)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        x = _eval_constexpr(node.operand, text, stack)
        return x if isinstance(node.op, ast.UAdd) else -x
    if isinstance(node, ast.BinOp):
        a = _eval_constexpr(node.left, text, stack)
        b = _eval_constexpr(node.right, text, stack)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            return a * b
        if isinstance(node.op, ast.Div):
            return a / b
    raise RuntimeError(f"unsupported constexpr expression node: {ast.dump(node)}")


def parse_const(text: str, name: str, stack: tuple[str, ...] = ()) -> float:
    """Resolve a scalar constexpr using only names, numbers and +,-,*,/."""
    if name in stack:
        raise RuntimeError(f"cyclic implementation constant alias: {' -> '.join((*stack, name))}")
    expressions = {n: expr for n, expr in CONST_RE.findall(text)}
    if name not in expressions:
        raise RuntimeError(f"cannot extract implementation constant {name}")
    expr = " ".join(_strip_float_suffixes(expressions[name]).split())
    try:
        tree = ast.parse(expr, mode="eval")
        value = _eval_constexpr(tree, text, (*stack, name))
    except (SyntaxError, ZeroDivisionError) as exc:
        raise RuntimeError(f"cannot evaluate implementation constant {name}: {expr!r}") from exc
    if not math.isfinite(value):
        raise RuntimeError(f"implementation constant {name} is non-finite: {expr!r}")
    return float(value)


def parse_aw_sigma_floor(text: str) -> float:
    """Extract the deployed lower floor used before setting Sigma_aw_stat."""
    pat = re.compile(
        r"const\s+float\s+sigma_floor\s*=\s*std::max\(\s*"
        r"([0-9.+\-eE]+)f?\s*,\s*band_noise_floor_sigma_\(\)\s*\)\s*;"
    )
    m = pat.search(text)
    if not m:
        raise RuntimeError("cannot extract deployed a_w stationary-std floor")
    return float(m.group(1))


def build(header: Path) -> dict:
    text = header.read_text()
    c = {name: parse_const(text, name) for name in REQUIRED}
    sigma_floor = parse_aw_sigma_floor(text)
    return {
        "schema": 2,
        "claim": "OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "source_complete_parameter_domain": True,
        "validated_arithmetic": False,
        "outward_rounded": False,
        "implementation_header": str(header.relative_to(REPO)),
        "continuous_parameters": {
            "wave_tune_frequency_hz": [c["MIN_TUNE_FREQ_HZ"], c["MAX_TUNE_FREQ_HZ"]],
            "tau_aw_s": [c["MIN_TAU_S"], c["MAX_TAU_S"]],
            "sigma_aw_mps2": [sigma_floor, c["MAX_SIGMA_A"]],
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
            "tilt_recovery": ["normal", "reset", "relock", "cooldown_reentry"],
            "aw_covariance_sync": ["not_due", "due_psd_increment"],
        },
        "hybrid_obligations": list(HYBRID_OBLIGATIONS),
        "periodic_aw_covariance_sync_proof": {
            "required_mode": "PSD_NONEXPANSIVE",
            "operation": "P_plus=P_minus+E_a Delta_plus E_a^T with Delta_plus>=0",
            "metric_consequence": "inverse-covariance information energy is nonexpansive",
        },
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
