#!/usr/bin/env python3
"""Fast refinement entry point for the physical-period r_S optimum sweep.

The first pass showed the vertical optimum was interior, but the 3-D optimum
reached the low-r_S edge.  Reuse the same experiment while compiling only the
simulator target instead of every OU-III unit-test binary.
"""

from __future__ import annotations

import subprocess

import ou_rs_adaptation_optimum_run as run


def _build_simulator_only(families, eigen_dir):
    if tuple(families) != ("OU_III",):
        raise ValueError("refinement runner is OU_III-only")
    command = ["make", "-B", "kalman_ou_iii-sim"]
    if eigen_dir:
        command.append(f"EIGEN_DIR={eigen_dir}")
    subprocess.run(
        command,
        cwd=run.core.FAMILY_MAKE_DIR["OU_III"],
        check=True,
    )


run.core.build_simulators = _build_simulator_only

if __name__ == "__main__":
    raise SystemExit(run.main())
