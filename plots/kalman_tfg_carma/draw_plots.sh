#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash draw_plots.sh`.
set -e

# The simulator writes its replays into tests/kalman_tfg_carma; the CI step
# copies them here before calling this. Running standalone, pick them up
# directly so the script is usable outside CI.
shopt -s nullglob
existing=( *_fusion_TFG_CARMA3_T41.csv )
if [ ${#existing[@]} -eq 0 ]; then
  produced=( ../../tests/kalman_tfg_carma/*_fusion_TFG_CARMA3_T41.csv )
  if [ ${#produced[@]} -gt 0 ]; then
    cp -f "${produced[@]}" .
  fi
fi

python3 kalman_tfg_carma-plots.py
