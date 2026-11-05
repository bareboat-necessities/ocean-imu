#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash draw_plots.sh`.
set -e

# The simulator writes its replays into tests/kalman_carma; the CI step copies
# them here before calling this. Running standalone, pick them up directly so
# the script is usable outside CI.
shopt -s nullglob
existing=( *_fusion_CARMA3_A41.csv )
if [ ${#existing[@]} -eq 0 ]; then
  produced=( ../../tests/kalman_carma/*_fusion_CARMA3_A41.csv )
  if [ ${#produced[@]} -gt 0 ]; then
    cp -f "${produced[@]}" .
  fi
fi

python3 kalman_carma-plots.py

# Coupled three-axis priors (doc/kalman_carma/kalman_carma_coupled.tex): the
# figures are drawn from the committed study data in
# reports/results/carma3_coupled (tools/carma3_coupled_benchmark.py).
python3 ../../tools/carma3_coupled_analysis.py --plots-from-json . --ext png
python3 ../../tools/carma3_coupled_report.py --figures-only --outdir . --ext png
