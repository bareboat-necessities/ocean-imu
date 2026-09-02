#!/usr/bin/env python3
"""Scaled-inverse augmented complete-word OU-III P4 design.

V4 reaches the augmented innovation solve but the raw-unit midpoint/Neumann
verifier can fail on a certified SPD 3x3 innovation family whose diagonal axes
have very different scales.  This wrapper adds an exact point diagonal
congruence preconditioner before giving up:

    T = D S D,              S^-1 = D T^-1 D,

where D is any positive point diagonal matrix.  D is chosen from binary64
midpoint diagonal magnitudes only for conditioning; no floating inverse is used
as an enclosure.  T is still an interval family exactly corresponding to S,
and its inverse is certified by the existing fixed-pivot interval
Gauss--Jordan / midpoint-Neumann verifier.  The recovered S^-1 is outward
interval arithmetic.

The augmented inverse derivative remains dS^-1=-S^-1(dS)S^-1 because the
existing augmented layer differentiates the returned verified enclosure after
this value solve.  No K interval matrix is materialized.

This remains a design gate and cannot promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_augmented_complete_word_design_v4 as V4

V3 = V4.V3
M = V4.M
JJ = V4.AUG.JJ
SCHEMA = 5
_ORIGINAL_VERIFIED_INVERSE = JJ.verified_inverse


def _diag_scale(S):
    d = []
    for i in range(len(S)):
        m = 0.5 * (float(S[i][i].lo) + float(S[i][i].hi))
        if not (math.isfinite(m) and m > 0.0):
            raise RuntimeError(f"innovation diagonal midpoint is not positive at {i}: {m}")
        # Any positive point scale gives an exact congruence identity.  The
        # sqrt choice is purely a numerical preconditioner.
        di = 1.0 / math.sqrt(m)
        if not (math.isfinite(di) and di > 0.0):
            raise RuntimeError(f"invalid innovation diagonal scale at {i}: {di}")
        d.append(di)
    return d


def scaled_verified_inverse(S):
    try:
        return _ORIGINAL_VERIFIED_INVERSE(S)
    except Exception as first:
        d = _diag_scale(S)
        n = len(S)
        T = [[JJ.I(d[i]) * S[i][j] * JJ.I(d[j]) for j in range(n)] for i in range(n)]
        try:
            Xt, meta = _ORIGINAL_VERIFIED_INVERSE(T)
        except Exception as second:
            raise type(second)(
                f"raw verified inverse failed ({type(first).__name__}: {first}); "
                f"diagonal-scaled verified inverse also failed: {second}"
            ) from second
        X = [[JJ.I(d[i]) * Xt[i][j] * JJ.I(d[j]) for j in range(n)] for i in range(n)]
        X = JJ.matrix_symmetric_hull(X)
        out = dict(meta)
        out.update({
            "inverse_backend": "SYMMETRIC_DIAGONAL_SCALED_" + str(meta.get("inverse_backend")),
            "diagonal_congruence_preconditioner_used": True,
            "diagonal_point_scales": d,
            "preconditioner_is_exact_point_congruence": True,
            "recovery_identity": "S^-1=D(DSD)^-1D",
            "ordinary_float_inverse_used_as_enclosure": False,
            "K_interval_matrix_materialized": False,
            "raw_inverse_failure": f"{type(first).__name__}: {first}",
        })
        return X, out


# AUG.verified_inverse_ad resolves JJ.verified_inverse at call time.
JJ.verified_inverse = scaled_verified_inverse


def build(domain_path: Path = M.DEFAULT_DOMAIN, *, source_node_index: int = 0,
          attitude_halfwidth_factor: float = 1.0 / 32.0,
          epsilon: float = 0.5):
    d = V4.build(
        domain_path,
        source_node_index=source_node_index,
        attitude_halfwidth_factor=attitude_halfwidth_factor,
        epsilon=epsilon,
    )
    d["schema_v5"] = SCHEMA
    d["qualification_v5"] = "OU3_P4_SCALED_VERIFIED_INVERSE_AUGMENTED_WORD_DESIGN"
    d["scaled_verified_inverse_available"] = True
    d["scaled_inverse_uses_only_point_congruence_and_verified_interval_inverse"] = True
    d["ordinary_float_inverse_used_as_enclosure"] = False
    d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"] = False
    d["P4_USABLE_CERTIFICATE_PROMOTED"] = False
    return d


def validate(d):
    f = V4.validate(d)
    if d.get("schema_v5") != SCHEMA:
        f.append("schema_v5 mismatch")
    for key in (
        "scaled_verified_inverse_available",
        "scaled_inverse_uses_only_point_congruence_and_verified_interval_inverse",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("ordinary_float_inverse_used_as_enclosure") is not False:
        f.append("ordinary floating inverse used as enclosure")
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=M.DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--attitude-halfwidth-factor", type=float, default=1.0 / 32.0)
    ap.add_argument("--epsilon", type=float, default=0.5)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(
        a.domain,
        source_node_index=a.source_node_index,
        attitude_halfwidth_factor=a.attitude_halfwidth_factor,
        epsilon=a.epsilon,
    )
    vf = validate(d)
    d["validation_pass_v5"] = not vf
    d["validation_failures_v5"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "beta_measurement": d.get("modes", {}).get(m, {}).get("beta_measurement_design_upper"),
            "operations": d.get("modes", {}).get(m, {}).get("operation_count"),
            "inverse_backends": [op.get("inverse_backend") for op in d.get("modes", {}).get(m, {}).get("operations",[])],
            "max_raw_correction": max([float(op.get("correction_theta_norm_upper_raw",0.0)) for op in d.get("modes",{}).get(m,{}).get("operations",[])] or [0.0]),
            "max_tight_correction": max([float(op.get("correction_theta_norm_upper",0.0)) for op in d.get("modes",{}).get(m,{}).get("operations",[])] or [0.0]),
        } for m in ("H","A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
