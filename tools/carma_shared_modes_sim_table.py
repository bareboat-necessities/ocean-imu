#!/usr/bin/env python3
"""Generate the deterministic last-900 s simulation table for the shared-modes QMEKF.

The shared-modes counterpart of ``tools/carma_sim_table.py``: the same eight
records in the same order as the OU-III, TFG and CARMA3 tables, the same
columns (X, Y and Z error as a percentage of the incident H_s, the 3D error
and the attitude) and the same deterministic protocol: one realization per
record, default seeds, the final 900 s of a 20-minute replay of the default
variant CSM.  The paired multi-seed comparison with CARMA3 is
reports/results/kalman_carma_shared_modes/.

Usage:
    python3 tools/carma_shared_modes_sim_table.py \\
        --binary tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim \\
        --data-dir plots/kalman_carma_shared_modes \\
        --output doc/kalman_carma_shared_modes/csm-sim-results-generated.tex-part
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from carma_sim_table import RECORDS, format_table, run_record  # noqa: E402

CSM_TABLE = {
    "tool": "tools/carma_shared_modes_sim_table.py",
    "variant": "shared-modes QMEKF (variant CSM)",
    "caption": "Shared-modes QMEKF (shared directional CARMA(4,1) modes with coupling confidence) "
               "simulation results (last \\SI{900}{s} RMS window).",
    "label": "tab:csm_sim_results",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim")
    parser.add_argument("--data-dir", default="plots/kalman_carma_shared_modes")
    parser.add_argument("--output", default="doc/kalman_carma_shared_modes/csm-sim-results-generated.tex-part")
    args = parser.parse_args()

    def resolve(p: str) -> Path:
        return Path(p) if Path(p).is_absolute() else REPO_ROOT / p

    binary, data_dir, output = resolve(args.binary), resolve(args.data_dir), resolve(args.output)
    if not binary.exists():
        print(f"simulator not found: {binary}", file=sys.stderr)
        return 1
    rows: list[tuple[str, float, dict[str, Any]]] = []
    for label, hs, name in RECORDS:
        record = data_dir / name
        if not record.exists():
            print(f"record not found: {record}", file=sys.stderr)
            return 1
        rows.append((label, hs, run_record(binary, record)))
        print(f"scored {name}", file=sys.stderr)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(format_table(rows, CSM_TABLE), encoding="utf-8")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
