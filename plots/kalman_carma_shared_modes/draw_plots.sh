#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash draw_plots.sh`.
set -e

# Figures of doc/kalman_carma_shared_modes, drawn from the committed
# validation data in reports/results/kalman_carma_shared_modes
# (tools/carma_shared_modes_validation.py) as PGF for the articles and SVG
# for review.  The dataset-replay tracking charts use the replay CSVs this
# directory receives from tests/kalman_carma_shared_modes when present and
# the committed traces otherwise.
python3 ../../tools/carma_shared_modes_validation.py figures --outdir . --ext pgf,svg
