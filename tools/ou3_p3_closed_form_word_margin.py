#!/usr/bin/env python3
"""Closed-form interval-certified whole-word translation margin.

Both sides of the P3 comparison are built in closed form from one observability
Gramian referenced to the covariance-word **endpoint**, so neither side steps a
Riccati recursion.  That matters: this ledger records three independent failures
of stepwise interval covariance propagation over a word, the last of which lost
positivity between 343 and 503 of 635 samples even in Joseph form and even with
the x cell narrowed 4096x.  Closed form is not an optimisation here, it is the
only route that survives interval arithmetic.

Ceiling (upper bound on P at the word endpoint)
-----------------------------------------------
Each S observation at lag ``ell`` from the endpoint constrains the endpoint state
through ``[ell^2/2, ell, 1]``, with its own process-corruption inflation
``R_eff = rmax + sigma^2 (ell^3/6)^2 + q_c ell^7/252``.  Then
``P_end <= (sum_k h_k h_k' / R_eff,k)^-1``.  No ``P0`` appears, so the ceiling is
independent of the initial covariance.

Floor (lower bound on P at the word endpoint)
---------------------------------------------
The floor must survive every *other* measurement channel.  The accelerometer
Jacobian has zero columns for ``v``, ``p`` and ``S`` -- ``measurement_update_acc_only``
forms its innovation covariance from ``OFF_TH``, ``OFF_AW``, ``OFF_BA`` and
``OFF_BG`` only -- so it reaches translation solely through the cross-covariance
with ``a_w``; the magnetometer informs attitude alone.  In the limiting case
where ``a_w`` is known *exactly* over the whole word, ``v(T) = v(0) + int a_w``
still carries ``v(0)``'s uncertainty, and likewise ``p`` and ``S``.  So the floor
credits the S observations at **full** strength (``R_eff = rmax``, no process
inflation, the most informative case and hence the smallest floor) and adds the
raw shipping prior information.  Propagation only injects process noise and so
only reduces information, giving ``Y0_propagated <= Y0`` and therefore
``(Y0 + G_full)^-1 <= truth``.

Within one source cell the two sides take opposite endpoints -- the ceiling the
fewest observations and weakest measurement, the floor the most and strongest --
so both hold for every parameter value the cell admits.

Non-promotion
-------------
This computes a margin.  It is **not** a P3 verdict: the canonical producer and
gate remain the promotion authority, and this module cannot set a theorem flag.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 1
WORD_SAMPLES = 635
GATE = 1.0e-18
CHANNELS = ("v", "p", "S")

# Shipping initial covariance, Kalman3D_Wave_OU_III constructor.
SIGMA_V0, SIGMA_P0, SIGMA_S0 = 1.0, 20.0, 50.0


def _I(x: float) -> Interval:
    return Interval.point(float(x))


def _consts() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")

    def c(name: str) -> float:
        m = re.search(r"constexpr float\s+" + name + r"\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"shipping constant {name} not found")
        return float(m.group(1))

    # Cadence constants come from the repo's own source parser rather than a
    # second regex here: PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT is FREQ_SMOOTHER_DT,
    # not a literal, and source_schedule() already resolves the clamps and the
    # constexpr-float ratio exactly as the shipping code computes them.
    sched = BASE.source_schedule()
    return {
        "rs_coeff": c("R_S_MSE_COEFF_DEFAULT"),
        "accel_density": c("R_S_ACCEL_NOISE_DENSITY_DEFAULT"),
        "min_rs": c("MIN_R_S"),
        "max_rs": c("MAX_R_S"),
        "period_min": float(sched["pseudo_min_s"]),
        "period_max": float(sched["pseudo_max_s"]),
        "ratio": float(sched["pseudo_ratio"]),
    }


def cadence_s(tau: float, k: dict) -> float:
    return min(max(k["ratio"] * float(tau), k["period_min"]), k["period_max"])


def rs_target(tau: float, sigma: float, k: dict) -> float:
    """Deployed SpectralMSE R_S, clamped exactly as the shipping code clamps it."""
    q = (2.0 * k["accel_density"]) ** (1.0 / 14.0)
    val = (k["rs_coeff"] * q * ((sigma * tau ** 4) ** (6.0 / 7.0))
           / math.sqrt(cadence_s(tau, k)))
    return min(max(val, k["min_rs"]), k["max_rs"])


def _inv3(G):
    """Interval 3x3 inverse; None when the determinant is not validated positive."""
    d = (G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1])
         - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
         + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]))
    if d.lo <= 0.0 <= d.hi:
        return None
    C = [[G[(i + 1) % 3][(j + 1) % 3] * G[(i + 2) % 3][(j + 2) % 3]
          - G[(i + 1) % 3][(j + 2) % 3] * G[(i + 2) % 3][(j + 1) % 3]
          for j in range(3)] for i in range(3)]
    di = d.reciprocal()
    return [[C[j][i] * di for j in range(3)] for i in range(3)]


def _gramian(tau: float, sigma: float, rs: float, k: dict, *, inflate: bool):
    """Endpoint-referenced S-observation Gramian over one covariance word."""
    h = float(BASE.source_schedule()["dt_s"])
    Tw = BASE.up(WORD_SAMPLES * h)
    ts = cadence_s(tau, k)
    n = int(Tw / ts)
    if n < 3:
        return None, n
    fac = max(BASE.source_schedule()["R_S_axis_std_factors"])
    rmax = BASE.up((rs * fac) ** 2)
    qc = BASE.up(2.0 * sigma * sigma / tau)
    G = [[Interval(0.0, 0.0) for _ in range(3)] for _ in range(3)]
    for j in range(1, n + 1):
        ell = Tw - j * ts
        reff = rmax
        if inflate:
            reff = BASE.up(rmax + sigma * sigma * (ell ** 3 / 6.0) ** 2
                           + qc * abs(ell) ** 7 / 252.0)
        li = _I(BASE.up(ell))
        hv = [li.square() * _I(0.5), li, _I(1.0)]
        ri = _I(reff).reciprocal()
        for a in range(3):
            for b in range(3):
                G[a][b] = G[a][b] + (hv[a] * hv[b]) * ri
    return G, n


def node_margin(node: dict, k: dict) -> dict:
    """Margin for one source node, taking opposite cell endpoints on each side."""
    tau_lo, tau_hi = float(node["tau_s"][0]), float(node["tau_s"][1])
    sig_lo, sig_hi = (float(node["sigma_filter_committed_mps2"][0]),
                      float(node["sigma_filter_committed_mps2"][1]))
    # Ceiling: fewest observations and weakest S measurement in the cell.
    rs_ceil = max(rs_target(tau_hi, sig_hi, k), rs_target(tau_hi, sig_lo, k))
    Gc, n_ceil = _gramian(tau_hi, sig_hi, rs_ceil, k, inflate=True)
    # Floor: most observations and strongest S measurement in the cell.
    rs_floor = min(rs_target(tau_lo, sig_lo, k), rs_target(tau_lo, sig_hi, k))
    Gf, n_floor = _gramian(tau_lo, sig_lo, rs_floor, k, inflate=False)
    if Gc is None or Gf is None:
        return {"source_node": int(node["index"]), "validated": False,
                "reason": "fewer than three S observations in the word"}
    for a, v0 in enumerate((SIGMA_V0 ** 2, SIGMA_P0 ** 2, SIGMA_S0 ** 2)):
        Gf[a][a] = Gf[a][a] + _I(1.0 / v0)
    Ci, Fi = _inv3(Gc), _inv3(Gf)
    if Ci is None or Fi is None:
        return {"source_node": int(node["index"]), "validated": False,
                "reason": "interval inversion not validated"}
    ratios = [abs(Fi[i][i].lo) / abs(Ci[i][i].hi) for i in range(3)]
    idx = ratios.index(min(ratios))
    return {"source_node": int(node["index"]), "validated": True,
            "tau_s": [tau_lo, tau_hi], "R_S_ceiling": rs_ceil, "R_S_floor": rs_floor,
            "S_observations_ceiling": n_ceil, "S_observations_floor": n_floor,
            "channel_ratios": ratios, "margin": ratios[idx],
            "binding_channel": CHANNELS[idx]}


def build(domain_path: Path = DEFAULT_DOMAIN, *, stride: int = 1) -> dict:
    k = _consts()
    nodes = NODES.build()["nodes"]
    rows = [node_margin(n, k) for n in nodes[::max(1, int(stride))]]
    ok = [r for r in rows if r.get("validated")]
    worst = min(ok, key=lambda r: r["margin"]) if ok else None
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_CLOSED_FORM_WHOLE_WORD_TRANSLATION_MARGIN",
        "non_promoting": True,
        "certifies_theorem_stage": False,
        "interval_certified": True,
        "stepwise_riccati_recursion_used": False,
        "initial_covariance_from_shipping_constructor": [SIGMA_V0, SIGMA_P0, SIGMA_S0],
        "R_S_from_deployed_clamped_SpectralMSE_law": True,
        "canonical_gate": GATE,
        "word_samples": WORD_SAMPLES,
        "nodes_evaluated": len(rows),
        "nodes_validated": len(ok),
        "stride": int(stride),
        "worst": worst,
        "worst_margin": None if worst is None else worst["margin"],
        "worst_clears_canonical_gate": bool(worst and worst["margin"] > GATE),
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema changed")
    if d.get("non_promoting") is not True:
        f.append("margin producer must remain non-promoting")
    if d.get("certifies_theorem_stage") is not False:
        f.append("margin producer must not claim a theorem stage")
    if d.get("stepwise_riccati_recursion_used") is not False:
        f.append("closed-form claim broken: a stepwise recursion was used")
    if float(d.get("canonical_gate", 0.0)) != GATE:
        f.append("canonical usefulness gate must remain exactly 1e-18")
    if int(d.get("word_samples", 0)) != WORD_SAMPLES:
        f.append("covariance word sample count changed")
    if d.get("nodes_validated", 0) != d.get("nodes_evaluated", -1):
        f.append("some source nodes did not validate")
    for r in d.get("rows") or []:
        if r.get("validated") and not (math.isfinite(r["margin"]) and r["margin"] > 0.0):
            f.append(f"node {r['source_node']} produced a non-positive margin")
            break
    if not d.get("rows"):
        f.append("no source nodes evaluated")
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, stride=a.stride)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({x: y for x, y in d.items() if x != "rows"}, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
