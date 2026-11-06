#!/usr/bin/env python3
"""Charts for the shared-modes QMEKF's deterministic dataset replays.

Reads the ``*_fusion_CARMA_SHARED_MODES_CSM.csv`` replays that
``tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim`` writes
(default variant CSM: shared directional modes with the coupling confidence)
and draws the CARMA3 gallery of ``plots/kalman_carma/kalman_carma-plots.py``
for them, so the CSM, CARMA3 and OU-III galleries compare page by page.  CI
moves the PGF/SVG files into ``doc/kalman_carma_shared_modes/`` for
``kalman_carma_shared_modes-charts.tex``.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

_CARMA = Path(__file__).resolve().parents[1] / "kalman_carma" / "kalman_carma-plots.py"
_spec = importlib.util.spec_from_file_location("kalman_carma_plots", _CARMA)
carma_plots = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(carma_plots)

FAMILY = "CARMA_SHARED_MODES_CSM"
PREFIX = "w3d_csm"

if __name__ == "__main__":
    carma_plots.main(family=FAMILY, prefix=PREFIX)
