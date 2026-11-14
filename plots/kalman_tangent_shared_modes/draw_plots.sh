#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Publication figures come from the pinned six-seed study, not the unrelated
# single-seed replay. Reproduce them with the commands in the study README.
python3 "${script_dir}/../../tools/tangent_shared_modes_report.py" --check-only --allow-pinned-sources --fetch-pinned-source
