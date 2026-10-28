#!/usr/bin/env python3
"""Charts for the TFG-CARMA3 filter's deterministic dataset replays.

Reads the ``*_fusion_TFG_CARMA3_T41.csv`` replays that
``tests/kalman_tfg_carma/tfg_carma3-sim`` writes (CARMA(4,1) on x, y and z on
the two-frame Lie group) and draws the CARMA3 gallery of
``plots/kalman_carma/kalman_carma-plots.py`` for them, so the two galleries
compare page by page.  CI moves the PGF/SVG files into ``doc/kalman_tfg_carma/``
for ``kalman_tfg_carma-charts.tex`` and ``kalman_tfg_carma.tex``.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_CARMA = Path(__file__).resolve().parents[1] / "kalman_carma" / "kalman_carma-plots.py"
_spec = importlib.util.spec_from_file_location("kalman_carma_plots", _CARMA)
carma_plots = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(carma_plots)

FAMILY = "TFG_CARMA3_T41"
PREFIX = "w3d_tfg_carma3"

if __name__ == "__main__":
    carma_plots.main(family=FAMILY, prefix=PREFIX)
