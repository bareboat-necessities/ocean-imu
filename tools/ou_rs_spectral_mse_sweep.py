#!/usr/bin/env python3
"""Sweep the bias-variance (SpectralMSE) integral-regularizer law.

The pole-placement derivation answers "what r_S preserves a chosen normalized
regularization pole kappa?".  It does not answer "what kappa minimizes
displacement error for a sea of amplitude sigma_a?", and that second question
is where the wave amplitude belongs: the pseudo-measurement both suppresses
low-frequency integration drift and distorts the genuine wave displacement, and
only the second cost scales with physical wave energy.

    J(omega_R) = 3 q_eff / (2 omega_R^3)  +  (1/2pi) int |G(jw)-1|^2 S_eta dw
    |G(jw)-1|^2 = (1 + 4x^2)/(1 + x^6),   x = w/omega_R

Well below the wave band |G-1|^2 -> 4/x^4, so the wave term is 4 m_-4 omega_R^4
and the balance gives omega_R*^7 = (9/32) q_eff / m_-4.  With m_-4 ~
sigma_a^2 tau^8 for a self-similar sea,

    r_S* = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S)

which is tau^(41/14) away from cadence clamps, against main's sigma_a tau^(5/2).

C_J absorbs (int x^-8 Phi_a(x) dx)^(3/7).  Evaluating the exact balance on the
eight reference spectra gives C_J ~ 0.054, so the centre of this sweep is an
analytical prediction rather than a fitted starting point; the sweep decides
whether it is also the complete-MEKF optimum.

The c_sigma axis is a control.  The law divides c_sigma back out to recover the
physical sigma_a,B, so r_S should be invariant to it and only the OU prior
should move.  A flat c_sigma axis confirms the two channels are cleanly
separated in this law.

Typical use:

    python3 tools/ou_rs_spectral_mse_sweep.py --jobs 4
"""

from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ou_rs_amplitude_retune_sweep import (  # noqa: E402
    DEFAULT_BINARY, DEFAULT_RECORD_DIR, RECORDS, run_record, summarize,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

# Analytical prediction from the exact balance on the eight reference spectra.
C_J_ANALYTICAL = 0.0538
C_J_GRID = (0.030, 0.038, 0.046, 0.0538, 0.062, 0.072, 0.085, 0.100)
C_SIGMA_GRID = (0.6, 0.9, 1.3, 1.8)


def parse_grid(text, default):
    if not text:
        return default
    return tuple(float(v) for v in text.replace(",", " ").split())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cj-grid")
    ap.add_argument("--c-sigma-grid")
    ap.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("--record-dir", type=Path, default=DEFAULT_RECORD_DIR)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_rs_spectral_mse")
    args = ap.parse_args(argv)

    cjs = parse_grid(args.cj_grid, C_J_GRID)
    csigs = parse_grid(args.c_sigma_grid, C_SIGMA_GRID)
    binary = args.binary.resolve()
    if not binary.exists():
        raise SystemExit(f"simulator not built: {binary}")
    record_dir = args.record_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_rows, grid_rows = [], []
    total = len(csigs) * len(cjs)
    done = 0
    for c_sigma in csigs:
        for c_J in cjs:
            env_extra = {"OU_III_RS_LAW": "3",
                         "OU_III_RS_MSE_COEFF": repr(c_J),
                         "OU_SIGMA_COEFF": repr(c_sigma)}
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                futures = [pool.submit(run_record, binary, record_dir / n, env_extra)
                           for _, _, n in RECORDS]
                results = [f.result() for f in futures]
            rows = [{"c_sigma": c_sigma, "c_J": c_J, "family": fam, "hs_m": hs,
                     "record": n, **m} for (fam, hs, n), m in zip(RECORDS, results)]
            all_rows.extend(rows)
            summary = {"c_sigma": c_sigma, "c_J": c_J, **summarize(rows)}
            grid_rows.append(summary)
            done += 1
            print(f"[{done:3d}/{total}] c_sigma={c_sigma:<5g} C_J={c_J:<7g} "
                  f"mean_z={summary['mean_z_pct_hs']:.4f} "
                  f"max_z={summary['max_z_pct_hs']:.4f} "
                  f"mean_3d={summary['mean_rms_3d_m']:.5f}", flush=True)

    for name, rows in (("raw.csv", all_rows), ("grid.csv", grid_rows)):
        with (args.output_dir / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    grid_rows.sort(key=lambda r: r["mean_z_pct_hs"])
    print(f"\nwrote {args.output_dir}/raw.csv and grid.csv\n")
    print("best 10 by mean vertical RMS (%Hs):")
    print(f"{'c_sigma':>8} {'C_J':>8} {'mean_z':>9} {'max_z':>9} {'mean_3d':>9}")
    for r in grid_rows[:10]:
        print(f"{r['c_sigma']:8g} {r['c_J']:8g} {r['mean_z_pct_hs']:9.4f} "
              f"{r['max_z_pct_hs']:9.4f} {r['mean_rms_3d_m']:9.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
