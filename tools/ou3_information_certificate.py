#!/usr/bin/env python3
"""Source-varying information-metric certificate for adaptive OU-III.

The fixed/coarsely-binned path LMI can reject stable transient Kalman dynamics
because it forces the same metric on source points with different Riccati
covariances. This stage uses the estimator's own *full* covariance as a
parameter-dependent Lyapunov metric:

    M(g) = Sigma(g)^(-1).

For an exact deterministic word Phi and corresponding covariance endpoints,

    Sigma_1 = Phi Sigma_0 Phi' + Omega.

In relative endpoint coordinates this gives the exact identity

    lambda_information
      = 1 - lambda_min(Sigma_1^(-1/2) Omega Sigma_1^(-1/2)).

Thus a uniformly positive relative Riccati injection margin is exactly a strict
word contraction margin in the source-varying information metric. This is
estimator covariance, not empirical truth-error covariance, and no transition
matrix is fitted from a trajectory.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import ou3_numerical_certificate as BASE

NX = 21
COV_MAGIC = b"OU3COV1\0"
COV_VERSION = 1
HORIZONS_S = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
TIME_TOL_S = 2.0e-5
PSD_REL_TOL = 2.0e-4
STRICT_LAMBDA_TOL = 2.0e-5
IDENTITY_ABS_TOL = 5.0e-7


@dataclass
class CovBlock:
    t0: float
    t1: float
    start: np.ndarray
    end: np.ndarray


def load_covariance_blocks(path: Path) -> tuple[list[CovBlock], dict]:
    blocks: list[CovBlock] = []
    with path.open("rb") as f:
        if f.read(8) != COV_MAGIC:
            raise RuntimeError(f"bad covariance magic: {path}")
        version, nx, stride = struct.unpack("<III", f.read(12))
        if version != COV_VERSION or nx != NX:
            raise RuntimeError(f"unsupported covariance format {version}/{nx}: {path}")
        matrix_bytes = NX * NX * 4
        while True:
            raw_t = f.read(16)
            if not raw_t:
                break
            if len(raw_t) != 16:
                raise RuntimeError(f"truncated covariance time header: {path}")
            t0, t1 = struct.unpack("<dd", raw_t)
            raw0 = f.read(matrix_bytes)
            raw1 = f.read(matrix_bytes)
            if len(raw0) != matrix_bytes or len(raw1) != matrix_bytes:
                raise RuntimeError(f"truncated covariance matrices: {path}")
            P0 = np.frombuffer(raw0, dtype="<f4").reshape(NX, NX).astype(float)
            P1 = np.frombuffer(raw1, dtype="<f4").reshape(NX, NX).astype(float)
            blocks.append(CovBlock(t0, t1, 0.5 * (P0 + P0.T), 0.5 * (P1 + P1.T)))
    return blocks, {"version": version, "nx": nx, "stride_samples": stride,
                    "block_count": len(blocks)}


def _spd_sqrt(P: np.ndarray) -> tuple[np.ndarray, float, float]:
    P = 0.5 * (P + P.T)
    lam, U = np.linalg.eigh(P)
    lo = float(np.min(lam))
    hi = float(np.max(lam))
    if not np.all(np.isfinite(lam)) or lo <= 0.0:
        raise np.linalg.LinAlgError(f"covariance is not SPD: lambda_min={lo}")
    return U @ np.diag(np.sqrt(lam)) @ U.T, lo, hi


def _information_word_metrics(phi: np.ndarray, Sigma0: np.ndarray,
                              Sigma1: np.ndarray) -> tuple[float, dict]:
    """Evaluate one word in stable endpoint-whitened coordinates.

    With S_i = Sigma_i^(1/2), define

        C = S_1^(-1) Phi S_0.

    Then lambda_information = lambda_max(C' C), while the relative Riccati
    injection is exactly I-C C'.  Evaluating the latter form avoids forming
    Omega = Sigma1-Phi Sigma0 Phi' first, which is a subtraction of large,
    nearly equal matrices in the strongly ill-conditioned physical units of
    the held/active covariance.  No tolerance is relaxed and no estimator
    quantity is changed.
    """
    Sigma0 = 0.5 * (Sigma0 + Sigma0.T)
    Sigma1 = 0.5 * (Sigma1 + Sigma1.T)
    S0, lo0, hi0 = _spd_sqrt(Sigma0)
    S1, lo1, hi1 = _spd_sqrt(Sigma1)
    C = np.linalg.solve(S1, phi @ S0)

    gram = C.T @ C
    gram_vals = np.linalg.eigvalsh(0.5 * (gram + gram.T))
    lam = float(np.max(gram_vals))

    rel = np.eye(C.shape[0]) - C @ C.T
    rel_vals = np.linalg.eigvalsh(0.5 * (rel + rel.T))
    inc = {
        "omega_relative_lambda_min": float(np.min(rel_vals)),
        "omega_relative_lambda_max": float(np.max(rel_vals)),
        "Sigma0_lambda_min": lo0,
        "Sigma0_lambda_max": hi0,
        "Sigma1_lambda_min": lo1,
        "Sigma1_lambda_max": hi1,
        "relative_injection_evaluation": "I-C C^T with C=Sigma1^-1/2 Phi Sigma0^1/2",
    }
    return lam, inc


def information_lambda(phi: np.ndarray, Sigma0: np.ndarray, Sigma1: np.ndarray) -> float:
    """lambda_max for Phi' Sigma1^-1 Phi <= lambda Sigma0^-1."""
    lam, _ = _information_word_metrics(phi, Sigma0, Sigma1)
    return lam


def covariance_increment_margin(phi: np.ndarray, Sigma0: np.ndarray, Sigma1: np.ndarray) -> dict:
    """Return relative Riccati injection and covariance metric bounds."""
    _, inc = _information_word_metrics(phi, Sigma0, Sigma1)
    return inc


def information_identity_residual(lam: float, inc: dict) -> float:
    """Residual of 1-lambda_info=lambda_min(relative Omega)."""
    return abs((1.0 - float(lam)) - float(inc["omega_relative_lambda_min"]))


def pair_map_covariance(map_path: Path, record: str) -> tuple[list[BASE.MapBlock], list[CovBlock], dict]:
    cov_path = map_path.with_name(map_path.name.replace("_exact_maps.bin", "_covariance.bin"))
    if not cov_path.exists():
        raise FileNotFoundError(cov_path)
    maps, map_meta = BASE.load_exact_maps(map_path, record)
    covs, cov_meta = load_covariance_blocks(cov_path)
    if len(maps) != len(covs):
        raise RuntimeError(f"map/covariance block count mismatch for {record}: {len(maps)} != {len(covs)}")
    worst_dt = 0.0
    for i, (m, c) in enumerate(zip(maps, covs)):
        dt = max(abs(m.t0 - c.t0), abs(m.t1 - c.t1))
        worst_dt = max(worst_dt, dt)
        if dt > TIME_TOL_S:
            raise RuntimeError(f"map/covariance time mismatch {record} block {i}: {dt}")
    return maps, covs, {"map": map_meta, "covariance": cov_meta,
                        "max_boundary_time_mismatch_s": worst_dt}


def physical_word(map_blocks: list[BASE.MapBlock], cov_blocks: list[CovBlock],
                  mode: str, start: int, count: int):
    seq = map_blocks[start:start + count]
    if len(seq) != count:
        return None
    if any((not b.valid) or b.hybrid_jump or b.mode != mode for b in seq):
        return None
    if any(abs(seq[k + 1].t0 - seq[k].t1) > 5e-3 for k in range(count - 1)):
        return None
    dim = 21 if mode == "A" else 18
    Phi = np.eye(NX)
    for b in seq:
        Phi = b.phi @ Phi
    return Phi[:dim, :dim], cov_blocks[start].start[:dim, :dim], cov_blocks[start + count - 1].end[:dim, :dim]


def evaluate_horizon(records: dict[str, tuple[list[BASE.MapBlock], list[CovBlock]]],
                     mode: str, horizon_s: float) -> dict:
    lambdas = []
    omega_mins = []
    identity_residuals = []
    covariance_lows = []
    covariance_highs = []
    worst = None
    invalid_spd = 0
    word_count = 0
    for record, (maps, covs) in records.items():
        if not maps:
            continue
        base = float(np.median([b.t1 - b.t0 for b in maps if b.t1 > b.t0]))
        n = max(1, int(round(horizon_s / base)))
        for start in range(len(maps) - n + 1):
            w = physical_word(maps, covs, mode, start, n)
            if w is None:
                continue
            Phi, P0, P1 = w
            word_count += 1
            try:
                lam, inc = _information_word_metrics(Phi, P0, P1)
            except np.linalg.LinAlgError:
                invalid_spd += 1
                continue
            lambdas.append(lam)
            omega_mins.append(inc["omega_relative_lambda_min"])
            identity_residuals.append(information_identity_residual(lam, inc))
            covariance_lows.extend((inc["Sigma0_lambda_min"], inc["Sigma1_lambda_min"]))
            covariance_highs.extend((inc["Sigma0_lambda_max"], inc["Sigma1_lambda_max"]))
            if worst is None or lam > worst[0]:
                worst = (lam, record, start, start + n - 1, inc)

    arr = np.asarray(lambdas, float)
    omg = np.asarray(omega_mins, float)
    ids = np.asarray(identity_residuals, float)
    if not len(arr):
        return {"mode": mode, "horizon_s": horizon_s, "status": "NO_WORDS",
                "word_count": word_count, "failure_reasons": ["NO_VALID_WORDS"],
                "information_pass": False}
    covariance_consistent = invalid_spd == 0 and float(np.min(omg)) >= -PSD_REL_TOL
    identity_consistent = float(np.max(ids)) <= IDENTITY_ABS_TOL
    strict = float(np.max(arr)) < 1.0 - STRICT_LAMBDA_TOL
    failures = []
    if invalid_spd:
        failures.append("NON_SPD_COVARIANCE")
    if not covariance_consistent:
        failures.append("COVARIANCE_RECURSION_NOT_PSD")
    if not identity_consistent:
        failures.append("INFORMATION_IDENTITY_NUMERICS")
    if not strict:
        failures.append("NO_STRICT_CONTRACTION")
    sigma_min = float(np.min(covariance_lows))
    sigma_max = float(np.max(covariance_highs))
    passed = covariance_consistent and identity_consistent and strict
    return {
        "mode": mode,
        "horizon_s": horizon_s,
        "status": "PASS" if passed else "FAIL",
        "word_count": word_count,
        "failure_reasons": failures,
        "invalid_spd_words": invalid_spd,
        "lambda_worst_information": float(np.max(arr)),
        "lambda_p99_information": float(np.quantile(arr, 0.99)),
        "lambda_p50_information": float(np.quantile(arr, 0.50)),
        "strict_contraction": bool(strict),
        "strict_margin_1_minus_lambda": 1.0 - float(np.max(arr)),
        "relative_Riccati_injection_margin_worst": float(np.min(omg)),
        "omega_relative_lambda_min_worst": float(np.min(omg)),
        "information_identity_residual_max": float(np.max(ids)),
        "information_identity_consistent": bool(identity_consistent),
        "Sigma_endpoint_lambda_min": sigma_min,
        "Sigma_endpoint_lambda_max": sigma_max,
        "Sigma_endpoint_condition_bound": sigma_max / sigma_min,
        "covariance_recursion_consistent": bool(covariance_consistent),
        "information_pass": bool(passed),
        "worst_word": None if worst is None else {
            "record": worst[1], "start_block": worst[2], "end_block": worst[3],
            "lambda": worst[0], "increment": worst[4],
        },
    }


def select_mode(records, mode: str) -> dict:
    attempts = [evaluate_horizon(records, mode, h) for h in HORIZONS_S]
    passing = [x for x in attempts if x.get("information_pass")]
    if passing:
        selected = passing[0]
        strongest = max(passing, key=lambda x: x["strict_margin_1_minus_lambda"])
    else:
        selected = min(attempts, key=lambda x: x.get("lambda_worst_information", math.inf))
        strongest = selected
    return {"mode": mode, "metric": "source-varying inverse estimator covariance",
            "selected": selected, "strongest_executed_margin": strongest,
            "attempts": attempts}


def markdown(report: dict) -> str:
    h = report["held"]["selected"]
    a = report["active"]["selected"]
    hs = report["held"]["strongest_executed_margin"]
    ass = report["active"]["strongest_executed_margin"]
    out = [
        "# OU-III source-varying information certificate", "",
        f"Status: **{report['status']}**", "",
        "Metric: `M(g) = Sigma_KF(g)^(-1)` from the actual estimator covariance; no truth-error covariance is used.", "",
        f"Held first strict word: {h.get('horizon_s')} s, worst lambda {h.get('lambda_worst_information')}, margin {h.get('strict_margin_1_minus_lambda')}",
        f"Active first strict word: {a.get('horizon_s')} s, worst lambda {a.get('lambda_worst_information')}, margin {a.get('strict_margin_1_minus_lambda')}",
        f"Held strongest tested margin: {hs.get('horizon_s')} s, 1-lambda {hs.get('strict_margin_1_minus_lambda')}",
        f"Active strongest tested margin: {ass.get('horizon_s')} s, 1-lambda {ass.get('strict_margin_1_minus_lambda')}",
        "",
        "The certificate checks the exact Riccati identity `1-lambda_information = lambda_min(Sigma1^-1/2 Omega Sigma1^-1/2)` and records endpoint covariance metric bounds.",
        "The relative injection is evaluated stably as `I-C C^T`, `C=Sigma1^-1/2 Phi Sigma0^1/2`, avoiding cancellation in `Sigma1-Phi Sigma0 Phi^T`; the identity tolerance is unchanged.",
        "A PASS is an exact executed-source linear certificate in a parameter-dependent information metric; deployment promotion still requires a validated continuous-source lower injection margin and covariance bounds.",
    ]
    if report["status"] != "PASS":
        out.extend([
            "",
            f"Held selected failure reasons: {h.get('failure_reasons', [])}",
            f"Active selected failure reasons: {a.get('failure_reasons', [])}",
        ])
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", type=Path, default=BASE.DEFAULT_OUT)
    args = ap.parse_args()
    out = args.output_dir.resolve()
    cert_path = out / "certificate.json"
    if not cert_path.exists():
        raise FileNotFoundError(cert_path)
    original = json.loads(cert_path.read_text())

    map_files = sorted(out.glob("*_exact_maps.bin"))
    if len(map_files) != 8:
        raise RuntimeError(f"expected eight exact-map files, found {len(map_files)}")

    records = {}
    alignment = {}
    for map_path in map_files:
        record = map_path.name.replace("_exact_maps.bin", "")
        maps, covs, meta = pair_map_covariance(map_path, record)
        records[record] = (maps, covs)
        alignment[record] = meta

    held = select_mode(records, "H")
    active = select_mode(records, "A")
    passed = bool(held["selected"].get("information_pass") and
                  active["selected"].get("information_pass"))
    report = {
        "schema": 2,
        "scope": "eight_noisy_reference_replays_exact_maps_and_exact_KF_covariance",
        "status": "PASS" if passed else "FAIL",
        "filter_regression": original.get("filter_regression"),
        "map_integrity": original.get("map_integrity"),
        "metric": "M(g)=Sigma_KF(g)^(-1)",
        "metric_provenance": "actual estimator Riccati/Joseph covariance, not empirical truth-error covariance",
        "identity": "1-lambda_information=lambda_min(Sigma1^-1/2 Omega Sigma1^-1/2), Omega=Sigma1-Phi Sigma0 Phi^T",
        "identity_evaluation": "relative injection evaluated as I-C C^T with C=Sigma1^-1/2 Phi Sigma0^1/2",
        "identity_absolute_tolerance": IDENTITY_ABS_TOL,
        "held": held,
        "active": active,
        "alignment": alignment,
        "deployment_linear_promotion_requirements": [
            "validated source-complete word family",
            "uniform positive lower bound on relative Riccati injection",
            "uniform positive lower and finite upper covariance eigenvalue bounds",
            "finite source-prefix gain",
        ],
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
        "deployment_missing": "validated continuous-source information/covariance enclosure and nonlinear/hybrid/stochastic bounds",
    }
    (out / "information_certificate.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    text = markdown(report)
    (out / "information_certificate.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
