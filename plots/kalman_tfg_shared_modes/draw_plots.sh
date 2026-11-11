#!/usr/bin/env bash
set -euo pipefail

# The report regenerates the bounded PGF/SVG figures and manuscript tables
# directly from the committed paired-study evidence. The article reads them
# from reports/results/tfg_shared_modes/figures; no dataset replay is needed.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python3 "${script_dir}/../../tools/tfg_shared_modes_report.py"
