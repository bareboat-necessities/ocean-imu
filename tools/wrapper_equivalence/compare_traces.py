#!/usr/bin/env python3
"""Compare two directories of wrapper characterization traces.

Every trace is a stream of raw float32 values written by the drivers in this
directory.  The comparison is exact: a trace passes only when it is
byte-identical to its reference.  For a trace that differs, the first
differing value and the largest absolute/relative difference are reported so
that a contraction-only drift (deployed flags) can be told apart from a
behavioral change.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def compare(ref: Path, new: Path) -> tuple[bool, str]:
    a = np.fromfile(ref, dtype="<f4")
    b = np.fromfile(new, dtype="<f4")
    if a.size != b.size:
        return False, f"length {a.size} != {b.size}"
    ua, ub = a.view("<u4"), b.view("<u4")
    diff = np.nonzero(ua != ub)[0]
    if diff.size == 0:
        return True, f"identical ({a.size} values)"
    finite = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a[finite].astype(np.float64) - b[finite].astype(np.float64))
    scale = np.maximum(np.abs(a[finite].astype(np.float64)), 1e-12)
    return False, (f"{diff.size} of {a.size} values differ, first at index {diff[0]} "
                   f"({a[diff[0]]!r} vs {b[diff[0]]!r}); max abs {d.max():.3e}, "
                   f"max rel {(d / scale).max():.3e}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", type=Path)
    ap.add_argument("candidate", type=Path)
    args = ap.parse_args()
    refs = sorted(args.reference.glob("*.bin"))
    if not refs:
        print("no reference traces", file=sys.stderr)
        return 2
    ok = True
    for ref in refs:
        new = args.candidate / ref.name
        if not new.is_file():
            print(f"MISSING {ref.name}")
            ok = False
            continue
        same, msg = compare(ref, new)
        ok &= same
        print(f"{'OK  ' if same else 'DIFF'} {ref.name}: {msg}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
