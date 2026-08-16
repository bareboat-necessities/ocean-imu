#!/usr/bin/env python3
"""Compatibility launcher for the full rS candidate validation experiment.

The original experiment runner formats an integral exponent as the invalid C++
literal ``3f``.  Keep the experiment definition unchanged and override only the
compile-macro formatter so integral-valued exponents are emitted as ``3.0f``.
"""

from __future__ import annotations

import os

import ou_rs_candidate_full_validation as exp


def build_arm(coeff: float | None, exponent: float | None) -> None:
    env = os.environ.copy()
    if coeff is None:
        env.pop("CPPFLAGS", None)
    else:
        exponent_literal = f"{float(exponent):.10g}"
        if "." not in exponent_literal and "e" not in exponent_literal.lower():
            exponent_literal += ".0"
        env["CPPFLAGS"] = " ".join(
            (
                "-DOU_III_EXPERIMENT_EFFECTIVE_RS_POWER=1",
                f"-DOU_III_EXPERIMENT_EFFECTIVE_RS_COEFF={coeff:.10g}f",
                f"-DOU_III_EXPERIMENT_EFFECTIVE_RS_TAU_EXP={exponent_literal}f",
            )
        )
    exp.run(
        (
            "make",
            "-C",
            str(exp.MAKE_DIR.relative_to(exp.REPO_ROOT)),
            "-B",
            "kalman_ou_iii-sim",
        ),
        env=env,
    )


exp.build_arm = build_arm

if __name__ == "__main__":
    raise SystemExit(exp.main())
