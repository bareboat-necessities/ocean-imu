#!/usr/bin/env python3
"""Rigorous promoted wrapper around the reusable PR #438 translation core.

PR #438 contains valuable complete-word exact-rational propagation, but its
``_mode`` setup formed several conditioned scale products in ordinary binary64
before wrapping them in intervals.  A rounded-up denominator can understate an
upper covariance bound, so that setup is not used for promotion here.

This wrapper keeps the exact-rational word propagation and 192-bit dyadic
Loewner compression from #438, while rebuilding every source conversion and
conditioned scale with outward interval arithmetic from the first operation.
C++ ``float`` literals are also rounded to their deployed binary32 values.

The result remains deliberately partial: it certifies complete-word
translation dissipation on the reconstructed limiting P3 source cell.  It does
not claim full-state nonlinear P4 or P5 completion.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import struct
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_translation_full_word_interval as CORE
import ou3_p4_translation_full_word_interval_fast as FAST  # patches CORE with dyadic Loewner compression

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = CORE.DEFAULT_DOMAIN
WRAPPER = CORE.WRAPPER
N = CORE.N


def _I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _source_float(text: str, name: str) -> float:
    """Return the exact deployed binary32 value of a C++ ``float`` literal."""
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot source-bind {name}")
    return struct.unpack("!f", struct.pack("!f", float(m.group(1))))[0]


def _mode(mode: str, domain_path: Path, horizon_s: float) -> dict:
    c = CORE.WORST.build_cell(mode, domain_path)
    s = CORE.WORST.serializable(c)
    row = c["row"]
    h = float(c["sched"]["dt_s"])
    x = c["x"]

    # No round-to-nearest division before intervalization.
    tau = _I(h) / x
    sigma = float(c["sigma"].lo)
    rho = float(c["rho_translation_lower"])

    # No round-to-nearest products before intervalization.  In particular the
    # covariance upper uses the *lower* possible scale^2 denominator through
    # interval division, so it cannot be understated by a rounded denominator.
    hs = _I(h)
    ss = _I(sigma)
    scales = [ss * hs, ss * hs * hs, ss * hs * hs * hs, ss]
    scale2 = [q.square() for q in scales]
    u = list(map(float, row["Sigma_diagonal_upper"]))
    physical = [u[6], u[9], u[12], u[15]]
    upper = [(_I(physical[i]) / scale2[i]).hi for i in range(N)]

    text = WRAPPER.read_text(encoding="utf-8")
    rh = min(
        _source_float(text, "R_S_x_factor_"),
        _source_float(text, "R_S_y_factor_"),
        1.0,
    )
    rs = (_I(rh) * _I(float(c["rs"].lo))).lo
    rsvar = _I(rs).square().lo
    acc = float(c["vector"]["configured_measurement_bounds"]["acc_measurement_std_mps2"])
    accvar = _I(acc).square().lo
    rS = (_I(rsvar) / scale2[2]).lo
    rA = (_I(accvar) / scale2[3]).lo
    betaS = (_I(1.0) / _I(rS)).hi
    betaA = (_I(1.0) / _I(rA)).hi

    n = int(math.ceil(horizon_s / h))
    leaves = CORE._prop(tau, h, n, rho, betaS, betaA)
    cert = []
    for t, L, radius, dep in leaves:
        d = CORE._delta(L, upper)
        if d <= 0.0 or not CORE._spd_delta(L, upper, d):
            raise RuntimeError(f"nonpositive endpoint margin on tau leaf {t.as_list()}")
        cert.append({
            "tau_s": t.as_list(),
            "delta_lower": d,
            "max_conditioned_radius_removed": radius,
            "split_depth": dep,
        })

    w = min(cert, key=lambda q: q["delta_lower"])
    old = float(row["direct_translation_generalized_margin_lower"])
    return {
        "source_cell": s,
        "conditioned_coordinates": "D^-1[v,p,S,a_w]",
        "tau_interval_s": tau.as_list(),
        "tau_leaf_count": len(cert),
        "max_tau_split_depth_used": max(q["split_depth"] for q in cert),
        "steps": n,
        "horizon_s": horizon_s,
        "process_injection_lower_conditioned": rho,
        "S_measurement_information_beta_conditioned": betaS,
        "accelerometer_aw_information_beta_conditioned": betaA,
        "measurement_information_geometry": "rank_one_S_and_aw_each_sample_exact_rational",
        "prediction_enclosure": "exact_rational_transition_interval_rowwise_loewner_dyadic192",
        "exact_rational_transition_enclosure": True,
        "dyadic_loewner_compression": True,
        "dyadic_loewner_bits": FAST.DYADIC_BITS,
        "source_float_literals_rounded_as_binary32": True,
        "conditioned_scale_products_outward_rounded_from_first_operation": True,
        "pre_interval_binary64_scale_products_used": False,
        "exact_rational_lower_retained_through_word": True,
        "corrections_allowed_every_sample_for_lower_bound": True,
        "artificial_S_variance_conditioned": rS,
        "artificial_acc_aw_variance_conditioned": rA,
        "translation_covariance_upper_conditioned": upper,
        "tau_leaf_certificates": cert,
        "complete_word_translation_margin_lower": w["delta_lower"],
        "limiting_tau_leaf": w,
        "old_single_seed_translation_margin_lower": old,
        "margin_widening_factor_lower": CORE.down(w["delta_lower"] / old),
        "interval_ldlt_endpoint_recertified": True,
    }


def build(domain_path=DEFAULT_DOMAIN, horizon_s=CORE.DEFAULT_HORIZON_S) -> dict:
    horizon_s = float(horizon_s)
    if not math.isfinite(horizon_s) or horizon_s <= 0.0:
        raise ValueError("horizon_s must be finite positive")
    p = Path(domain_path).resolve()
    modes, failures = {}, []
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode(mode, p, horizon_s)
        except Exception as exc:
            failures.append(f"{mode}: {exc}")
    return {
        "qualification": "OU3_P4_RIGOROUS_REUSED_438_COMPLETE_WORD_TRANSLATION",
        "reused_from_PR_438": True,
        "source_only": True,
        "trajectory_replay_used": False,
        "outward_rounded": True,
        "unsafe_438_pre_interval_scale_setup_used": False,
        "horizon_s": horizon_s,
        "modes": modes,
        "P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": (
            "PASS" if not failures and len(modes) == 2 else "NOT_ESTABLISHED"
        ),
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "remaining_obligation": (
            "compose this complete-word translation margin with operation-matched finite-angle attitude/bias dissipation on every source-reachable path cell"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("reused_from_PR_438") is not True:
        failures.append("PR #438 provenance missing")
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        failures.append("source qualification invalid")
    if d.get("outward_rounded") is not True:
        failures.append("promoted wrapper is not outward rounded")
    if d.get("unsafe_438_pre_interval_scale_setup_used") is not False:
        failures.append("unsafe #438 scale setup was reused")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if not float(m.get("complete_word_translation_margin_lower", 0.0)) > 0.0:
            failures.append(f"{mode}: no complete-word translation margin")
        if not float(m.get("margin_widening_factor_lower", 0.0)) > 1.0:
            failures.append(f"{mode}: full-word translation did not widen P3 seed")
        if m.get("source_float_literals_rounded_as_binary32") is not True:
            failures.append(f"{mode}: source float semantics not binary32")
        if m.get("conditioned_scale_products_outward_rounded_from_first_operation") is not True:
            failures.append(f"{mode}: scale products not outward rounded from first op")
        if m.get("pre_interval_binary64_scale_products_used") is not False:
            failures.append(f"{mode}: binary64 pre-interval scale products returned")
        if m.get("dyadic_loewner_compression") is not True:
            failures.append(f"{mode}: runtime-bounded exact Loewner compression missing")
        if m.get("interval_ldlt_endpoint_recertified") is not True:
            failures.append(f"{mode}: endpoint not recertified")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        failures.append("translation-only result prematurely promoted P4")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--horizon-s", type=float, default=CORE.DEFAULT_HORIZON_S)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain, args.horizon_s)
    failures = validate(out)
    out["validation_failures"] = failures
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS"],
        "modes": {m: {
            "delta": out.get("modes", {}).get(m, {}).get("complete_word_translation_margin_lower"),
            "factor": out.get("modes", {}).get(m, {}).get("margin_widening_factor_lower"),
            "tau_leaves": out.get("modes", {}).get(m, {}).get("tau_leaf_count"),
        } for m in ("H", "A")},
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())