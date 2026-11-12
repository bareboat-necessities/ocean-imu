#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"
if [[ $# -eq 1 && $1 == --study-only ]]; then
  # The comparative study uses committed scores and compact traces.
  python3 "${script_dir}/../../tools/tfg_shared_modes_report.py"
  exit 0
elif [[ $# -ne 0 ]]; then
  echo "Usage: $0 [--study-only]" >&2
  exit 2
fi
shopt -s nullglob
produced=( ../../tests/kalman_tfg_shared_modes/*_fusion_TFG_CARMA_SHARED_MODES_CSM.csv )
if [ ${#produced[@]} -gt 0 ]; then
  # Refresh earlier or partial plot inputs after every new replay.
  cp -f "${produced[@]}" .
fi
if ! python3 -c 'import glob; assert len(glob.glob("*_fusion_TFG_CARMA_SHARED_MODES_CSM.csv")) >= 8'; then
  echo "Run tests/kalman_tfg_shared_modes/run_tests.sh before drawing the replay gallery." >&2
  exit 1
fi
python3 kalman_tfg_shared_modes-plots.py
# Both the method article and the comparative study reuse these bounded
# figures and tables from the committed paired-study evidence.
python3 "${script_dir}/../../tools/tfg_shared_modes_report.py"
