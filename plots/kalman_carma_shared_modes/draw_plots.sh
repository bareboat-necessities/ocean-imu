#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash draw_plots.sh`.
set -e

# 1. The per-record gallery of the deterministic dataset replays
#    (kalman_carma_shared_modes-charts.tex): CI copies the replay CSVs here;
#    a local run takes them from the test directory.
shopt -s nullglob
existing=( *_fusion_CARMA_SHARED_MODES_CSM.csv )
if [ ${#existing[@]} -eq 0 ]; then
  produced=( ../../tests/kalman_carma_shared_modes/*_fusion_CARMA_SHARED_MODES_*.csv )
  if [ ${#produced[@]} -gt 0 ]; then
    cp -f "${produced[@]}" .
  fi
fi
python3 kalman_carma_shared_modes-plots.py

# 2. The validation figures of kalman_carma_shared_modes-validation.tex, drawn
#    from the committed validation data in reports/results/kalman_carma_shared_modes
#    (tools/carma_shared_modes_validation.py) and, for the tracking charts,
#    from the replays above.
python3 ../../tools/carma_shared_modes_validation.py figures --outdir . --ext pgf,svg

# 3. The figures of kalman_carma_shared_modes-minimal.tex (the 27-state
#    minimal realization), from the committed study data in
#    reports/results/kalman_carma_shared_modes_minimal (analysis.json from
#    tools/carma_shared_modes_minimal_analysis.py, summary.json and traces from
#    tools/carma_shared_modes_minimal_study.py).
python3 ../../tools/carma_shared_modes_minimal_analysis.py --figures-only --outdir . --ext pgf,svg
python3 ../../tools/carma_shared_modes_minimal_study.py figures --outdir . --ext pgf,svg
