#!/usr/bin/env python3
"""Generate the deterministic last-900 s simulation table for TFG-CARMA3.

The two-frame Lie-group counterpart of ``tools/carma_sim_table.py``: the same
eight records in the same order as the OU-III, TFG and CARMA3 tables, the same
columns (X, Y and Z error as a percentage of the incident H_s, the 3D error
and the attitude) and the same historical deterministic protocol: one
realization per record, default seeds, the final 900 s of a 20-minute replay
of the default variant T41.  The paired multi-seed comparison with CARMA3 and
OU-III is reports/results/tfg_carma3/.

Usage:

    python3 tools/tfg_carma_sim_table.py \\
        --binary tests/kalman_tfg_carma/tfg_carma3-sim \\
        --data-dir plots/kalman_tfg_carma \\
        --output doc/kalman_tfg_carma/tfg-carma-sim-results-generated.tex-part
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from carma_sim_table import RECORDS, format_table, run_record  # noqa: E402

TFG_CARMA3_TABLE = {
    "tool": "tools/tfg_carma_sim_table.py",
    "variant": "TFG-CARMA3 (variant T41)",
    "caption": "TFG-CARMA3 (CARMA(4,1) on the two-frame Lie group) simulation results "
               "(last \\SI{900}{s} RMS window).",
    "label": "tab:tfg_carma_sim_results",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="tests/kalman_tfg_carma/tfg_carma3-sim")
    parser.add_argument("--data-dir", default="plots/kalman_tfg_carma")
    parser.add_argument("--output",
                        default="doc/kalman_tfg_carma/tfg-carma-sim-results-generated.tex-part")
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
    output.write_text(format_table(rows, TFG_CARMA3_TABLE), encoding="utf-8")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
