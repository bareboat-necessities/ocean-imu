#!/usr/bin/env python3
"""Extract the OU-III continuous-source proof domain from implementation guards.

This producer intentionally does not infer bounds from the eight reference
trajectories. It parses the shipping implementation's safety/clamp constants
and records every discrete branch and hybrid transition the validated backend
must cover.

The implementation constants are C++ ``float`` constexpr expressions. For a
proof-domain boundary it is not sufficient to evaluate those expressions as
Python binary64. This parser therefore evaluates every literal and arithmetic
operation as IEEE-754 binary32, using exact rationals between operations and an
explicit nearest/ties-to-even rounding step. The returned Python float is an
exact binary64 representation of the deployed binary32 value.
"""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import json
import math
import re
import struct
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
_FLOAT32_MAX_BITS = 0x7F7FFFFF


def _strip_float_suffixes(expr: str) -> str:
    return re.sub(r"(?<=\d)[fF]\b", "", expr)


def _bits_to_positive_fraction(bits: int) -> Fraction:
    """Return the exact nonnegative finite binary32 value for ``bits``."""
    if not 0 <= bits <= _FLOAT32_MAX_BITS:
        raise ValueError(f"not a finite positive binary32 pattern: 0x{bits:08x}")
    exponent = (bits >> 23) & 0xFF
    mantissa = bits & 0x7FFFFF
    if exponent == 0:
        return Fraction(mantissa, 1 << 149)
    significand = (1 << 23) | mantissa
    power = exponent - 127 - 23
    if power >= 0:
        return Fraction(significand << power, 1)
    return Fraction(significand, 1 << (-power))


def _round_fraction_binary32(value: Fraction) -> Fraction:
    """Round an exact rational to finite binary32, nearest/ties-to-even.

    Python binary64 is used only to locate a nearby candidate. The final choice
    is made by exact rational distances over that candidate and its neighbours,
    so possible binary64 double-rounding cannot change the selected binary32.
    """
    if value == 0:
        return Fraction(0, 1)
    sign = -1 if value < 0 else 1
    x = abs(value)
    try:
        candidate = struct.unpack(">I", struct.pack(">f", float(x)))[0]
    except OverflowError as exc:
        raise RuntimeError(f"binary32 constant overflow for {value}") from exc
    candidate &= 0x7FFFFFFF
    if candidate > _FLOAT32_MAX_BITS:
        raise RuntimeError(f"binary32 constant overflow for {value}")

    choices: list[tuple[Fraction, int]] = []
    for bits in (candidate - 1, candidate, candidate + 1):
        if 0 <= bits <= _FLOAT32_MAX_BITS:
            choices.append((_bits_to_positive_fraction(bits), bits))
    if not choices:
        raise RuntimeError(f"cannot round binary32 constant {value}")

    def rank(item: tuple[Fraction, int]) -> tuple[Fraction, int]:
        exact, bits = item
        # Exact midpoint: ties-to-even means an even low significand bit.
        return (abs(exact - x), bits & 1)

    exact, _ = min(choices, key=rank)
    return exact if sign > 0 else -exact


def _literal_fraction(expr: str, node: ast.Constant) -> Fraction:
    token = ast.get_source_segment(expr, node)
    if token is None:
        token = repr(node.value)
    try:
        return Fraction(token)
    except (ValueError, ZeroDivisionError) as exc:
        raise RuntimeError(f"cannot parse constexpr literal {token!r}") from exc


def _eval_constexpr32(node: ast.AST, expr: str, text: str,
                      stack: tuple[str, ...]) -> Fraction:
    if isinstance(node, ast.Expression):
        return _eval_constexpr32(node.body, expr, text, stack)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return _round_fraction_binary32(_literal_fraction(expr, node))
    if isinstance(node, ast.Name):
        return parse_const_fraction(text, node.id, stack)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        x = _eval_constexpr32(node.operand, expr, text, stack)
        return x if isinstance(node.op, ast.UAdd) else -x
    if isinstance(node, ast.BinOp):
        a = _eval_constexpr32(node.left, expr, text, stack)
        b = _eval_constexpr32(node.right, expr, text, stack)
        if isinstance(node.op, ast.Add):
            exact = a + b
        elif isinstance(node.op, ast.Sub):
            exact = a - b
        elif isinstance(node.op, ast.Mult):
            exact = a * b
        elif isinstance(node.op, ast.Div):
            if b == 0:
                raise ZeroDivisionError("constexpr division by zero")
            exact = a / b
        else:
            raise RuntimeError(f"unsupported constexpr operator: {ast.dump(node.op)}")
        return _round_fraction_binary32(exact)
    raise RuntimeError(f"unsupported constexpr expression node: {ast.dump(node)}")


def parse_const_fraction(text: str, name: str,
                         stack: tuple[str, ...] = ()) -> Fraction:
    """Resolve a scalar ``constexpr float`` with deployed binary32 semantics."""
    if name in stack:
        raise RuntimeError(f"cyclic implementation constant alias: {' -> '.join((*stack, name))}")
    expressions = {n: expr for n, expr in CONST_RE.findall(text)}
    if name not in expressions:
        raise RuntimeError(f"cannot extract implementation constant {name}")
    expr = " ".join(_strip_float_suffixes(expressions[name]).split())
    try:
        tree = ast.parse(expr, mode="eval")
        return _eval_constexpr32(tree, expr, text, (*stack, name))
    except (SyntaxError, ZeroDivisionError) as exc:
        raise RuntimeError(f"cannot evaluate implementation constant {name}: {expr!r}") from exc


def parse_const(text: str, name: str, stack: tuple[str, ...] = ()) -> float:
    value = float(parse_const_fraction(text, name, stack))
    if not math.isfinite(value):
        raise RuntimeError(f"implementation constant {name} is non-finite")
    return value


def parse_aw_sigma_floor(text: str) -> float:
    """Extract the deployed lower floor used before setting Sigma_aw_stat."""
    pat = re.compile(
        r"const\s+float\s+sigma_floor\s*=\s*std::max\(\s*"
        r"([0-9.+\-eE]+)f?\s*,\s*band_noise_floor_sigma_\(\)\s*\)\s*;"
    )
    m = pat.search(text)
    if not m:
        raise RuntimeError("cannot extract deployed a_w stationary-std floor")
    return float(_round_fraction_binary32(Fraction(m.group(1))))


def _outward_point(x: float) -> list[float]:
    """Return an explicit binary64 enclosure of an exact deployed binary32."""
    x = float(x)
    if not math.isfinite(x):
        raise RuntimeError(f"non-finite source-domain endpoint {x!r}")
    return [math.nextafter(x, -math.inf), math.nextafter(x, math.inf)]


def _outward_box(lo: float, hi: float) -> list[float]:
    """Return an outward-rounded binary64 box containing two source endpoints."""
    lo = float(lo)
    hi = float(hi)
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
        raise RuntimeError(f"invalid source interval [{lo!r}, {hi!r}]")
    return [math.nextafter(lo, -math.inf), math.nextafter(hi, math.inf)]


def build(header: Path) -> dict:
    text = header.read_text()
    c = {name: parse_const(text, name) for name in REQUIRED}
    sigma_floor = parse_aw_sigma_floor(text)
    continuous = {
        "wave_tune_frequency_hz": [c["MIN_TUNE_FREQ_HZ"], c["MAX_TUNE_FREQ_HZ"]],
        "tau_aw_s": [c["MIN_TAU_S"], c["MAX_TAU_S"]],
        "sigma_aw_mps2": [sigma_floor, c["MAX_SIGMA_A"]],
        "R_S_base": [c["MIN_R_S"], c["MAX_R_S"]],
        "pseudo_update_period_s": [
            c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"], c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"],
        ],
    }
    timing = {
        "mag_delay": c["MAG_DELAY_SEC"],
        "online_tune_warmup": c["ONLINE_TUNE_WARMUP_SEC"],
    }
    parameter_box = {
        "qualification": "SOURCE_DERIVED_OUTWARD_ROUNDED_PARAMETER_BOX",
        "validated_arithmetic": True,
        "outward_rounded": True,
        "arithmetic_backend": "EXACT_BINARY32_SOURCE_PLUS_BINARY64_NEXTAFTER_OUTWARD",
        "continuous_parameters": {
            name: _outward_box(bounds[0], bounds[1])
            for name, bounds in continuous.items()
        },
        "timing_constants_s": {
            name: _outward_point(value)
            for name, value in timing.items()
        },
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
    }
    return {
        "schema": 2,
        "claim": "OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "source_complete_parameter_domain": True,
        # These top-level flags remain false until the whole H/A word is
        # propagated with validated matrix/transcendental arithmetic.
        "validated_arithmetic": False,
        "outward_rounded": False,
        "implementation_header": str(header.relative_to(REPO)),
        "implementation_scalar_semantics": {
            "type": "IEEE754_BINARY32",
            "rounding": "ROUND_TO_NEAREST_TIES_TO_EVEN_EACH_OPERATION",
            "evaluation": "EXACT_RATIONAL_THEN_BINARY32_ROUND",
        },
        "continuous_parameters": continuous,
        "timing_constants_s": timing,
        "validated_parameter_box": parameter_box,
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
            "the validated_parameter_box closes only the source-boundary arithmetic; theorem "
            "promotion still requires outward-rounded interval/Taylor-model propagation of "
            "the complete H/A Riccati words, nonlinear SO(3) remainder and remaining jumps"
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
