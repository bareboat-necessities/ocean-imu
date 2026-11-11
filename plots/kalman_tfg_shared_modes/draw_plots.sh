#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"
shopt -s nullglob
existing=( *_fusion_TFG_CARMA_SHARED_MODES_CSM.csv )
if [ ${#existing[@]} -eq 0 ]; then
  produced=( ../../tests/kalman_tfg_shared_modes/*_fusion_TFG_CARMA_SHARED_MODES_CSM.csv )
  if [ ${#produced[@]} -gt 0 ]; then
    cp -f "${produced[@]}" .
  fi
fi
if ! python3 -c 'import glob; assert len(glob.glob("*_fusion_TFG_CARMA_SHARED_MODES_CSM.csv")) >= 8'; then
  echo "Run tests/kalman_tfg_shared_modes/run_tests.sh before drawing the replay gallery." >&2
  exit 1
fi
python3 kalman_tfg_shared_modes-plots.py
# Both the method article and the comparative study reuse these bounded
# figures and tables from the committed paired-study evidence.
python3 "${script_dir}/../../tools/tfg_shared_modes_report.py"
