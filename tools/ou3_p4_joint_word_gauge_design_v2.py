#!/usr/bin/env python3
"""Corrected zero-rate source realization for the active OU-III P4 route.

The original co-gauged design used the wrong sign for the zero-rate
attitude/gyro-bias process cross covariance.  The source has B(t)=t I at
omega=0 and integral_0^h B(s) ds = h^2/2 I, so the q_b cross term is positive.
This module keeps only that source correction.  The current augmented
complete-word backend implements the exact S-selector cancellation directly,
so the retired process-global S-update monkey patch is no longer imported here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_p4_joint_word_gauge_design as G

DEFAULT_DOMAIN = G.DEFAULT_DOMAIN
_ORIGINAL = G._zero_rate_transition_process


def corrected_zero_rate_transition_process(mode: str, src: dict, domain: dict):
    F, Q, meta = _ORIGINAL(mode, src, domain)
    h = float(src["dt_s"])
    qb = float(PROCESS.build()["source_constants"]["gyro_bias_rw_variance_density"])
    qtb = qb * h * h / 2.0
    for ax in range(3):
        z = Interval.outward_bounds(qtb, qtb)
        Q[ax][3 + ax] = z
        Q[3 + ax][ax] = z
    Q = G.H._psd_tighten(Q)
    return F, Q, meta


# Install only inside processes that explicitly import this active source fix.
G._zero_rate_transition_process = corrected_zero_rate_transition_process


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = G.build(a.domain, a.source_node_index)
    d["zero_rate_attitude_bias_process_cross_sign"] = "POSITIVE_SOURCE_EXACT"
    d["legacy_wrong_process_cross_sign_used"] = False
    d["physical_S_true_forcing_must_enter_affine_b"] = True
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "ops": d.get("modes", {}).get(m, {}).get("operation_count"),
            "last_op": (d.get("modes", {}).get(m, {}).get("operations") or [None])[-1],
        } for m in ("H", "A")},
        "failures": d["failures"],
    }, indent=2, sort_keys=True))
    return 0 if not d["failures"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
