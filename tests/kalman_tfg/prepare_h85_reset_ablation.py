#!/usr/bin/env python3
"""Prepare one reset-only TFG Hs=8.5 ablation in a GitHub Actions workspace.

The committed production implementation remains untouched.  Each CI arm checks
out the same corrected branch, then this helper changes only the covariance
transport applied after the same finite nominal injection.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
KALMAN = ROOT / "src/kalman_tfg/Kalman3D_Wave_TFG.h"

BASE = """        MatrixNX Jreset;\n        build_reset_jacobian(correction, Jreset);\n        inject(Tangent(-correction));\n        P_ = (Jreset * Pj * Jreset.transpose()).eval();\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {n}")
    return text.replace(old, new, 1)


def patch_diagonal_only(text: str) -> str:
    new = """        MatrixNX Jreset;\n        build_reset_jacobian(correction, Jreset);\n        // Ablation: retain the per-3D-block right-Jacobian terms but remove\n        // every semidirect-product vector<-attitude cross block.\n        for (int off = OFF_BG; off < NX; off += 3)\n            Jreset.template block<3,3>(off, OFF_PHI).setZero();\n        inject(Tangent(-correction));\n        P_ = (Jreset * Pj * Jreset.transpose()).eval();\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""
    return replace_once(text, BASE, new, "diagonal-only reset ablation")


def patch_inverse(text: str) -> str:
    new = """        MatrixNX Jreset;\n        build_reset_jacobian(correction, Jreset);\n        // Ablation: transport with the inverse of the currently implemented\n        // reset Jacobian while keeping the exact same nominal injection.\n        const MatrixNX Jtransport = Jreset.inverse();\n        inject(Tangent(-correction));\n        P_ = (Jtransport * Pj * Jtransport.transpose()).eval();\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""
    return replace_once(text, BASE, new, "inverse reset ablation")


def patch_opposite_sign(text: str) -> str:
    new = """        MatrixNX Jreset;\n        // Ablation: evaluate the same full group Jacobian at -a rather than\n        // +a.  Nominal injection remains -a in every arm.\n        build_reset_jacobian(Tangent(-correction), Jreset);\n        inject(Tangent(-correction));\n        P_ = (Jreset * Pj * Jreset.transpose()).eval();\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""
    return replace_once(text, BASE, new, "opposite-sign reset ablation")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: prepare_h85_reset_ablation.py full|diagonal_only|inverse|opposite_sign",
              file=sys.stderr)
        return 2

    variant = sys.argv[1]
    allowed = {"full", "diagonal_only", "inverse", "opposite_sign"}
    if variant not in allowed:
        print(f"unknown variant {variant!r}", file=sys.stderr)
        return 2

    text = KALMAN.read_text()
    if variant == "diagonal_only":
        text = patch_diagonal_only(text)
    elif variant == "inverse":
        text = patch_inverse(text)
    elif variant == "opposite_sign":
        text = patch_opposite_sign(text)

    KALMAN.write_text(text)
    print(f"prepared TFG H8.5 reset ablation variant: {variant}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
