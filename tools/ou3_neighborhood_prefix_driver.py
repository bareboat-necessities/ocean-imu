#!/usr/bin/env python3
"""Run the OU-III nonlinear neighborhood diagnostic on exact source prefixes.

The nonlinear observer must replay from t=0 so the adaptive estimator reaches the
same source state at injection, but it does not need the unused tail of a
20-minute source record after the certified word endpoint.  This wrapper creates
deterministic per-record prefixes, then invokes ``ou3_neighborhood_diagnostic``
unchanged on those prefixes.

This is a runtime optimization only.  It changes neither estimator behavior nor
sampled-certificate qualification.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import ou3_numerical_certificate as BASE

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "tools" / "ou3_neighborhood_diagnostic.py"


def _first_float(line: str) -> float | None:
    token = line.split(",", 1)[0].strip()
    try:
        value = float(token)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def write_prefix(source: Path, destination: Path, stop_s: float) -> dict:
    """Copy complete CSV rows through the first timestamp strictly after stop_s."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    copied_rows = 0
    last_t = None
    crossed = False
    with source.open("r", encoding="utf-8", errors="strict") as src, destination.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        for line in src:
            dst.write(line)
            t = _first_float(line)
            if t is None:
                continue
            copied_rows += 1
            last_t = t
            if t > stop_s:
                crossed = True
                break
    if copied_rows == 0:
        raise RuntimeError(f"no numeric source rows found in {source}")
    if not crossed:
        raise RuntimeError(
            f"source {source} ended at {last_t} before required prefix {stop_s}s"
        )
    return {
        "source": source.name,
        "stop_s": float(stop_s),
        "last_copied_time_s": float(last_t),
        "numeric_rows": copied_rows,
    }


def replace_option(argv: list[str], option: str, value: str) -> list[str]:
    out: list[str] = []
    i = 0
    replaced = False
    while i < len(argv):
        if argv[i] == option:
            if i + 1 >= len(argv):
                raise ValueError(f"{option} requires a value")
            out.extend((option, value))
            i += 2
            replaced = True
        else:
            out.append(argv[i])
            i += 1
    if not replaced:
        out.extend((option, value))
    return out


def main() -> int:
    probe = argparse.ArgumentParser(add_help=False)
    probe.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    probe.add_argument("--data-dir", type=Path, default=BASE.DEFAULT_DATA_DIR)
    probe.add_argument("--diagnostic-dir", type=Path, required=True)
    probe.add_argument("--modes", default="H,A")
    probe.add_argument("--held-time-s", type=float, default=60.0)
    probe.add_argument("--active-time-s", type=float, default=300.0)
    args, _ = probe.parse_known_args()

    modes = [x.strip() for x in args.modes.split(",") if x.strip()]
    if not modes or any(x not in ("H", "A") for x in modes):
        raise ValueError("--modes accepts H,A")
    contract = json.loads(
        (args.certificate_dir.resolve() / "information_enclosure_contract.json").read_text()
    )
    required_end = []
    for mode in modes:
        inject = args.held_time_s if mode == "H" else args.active_time_s
        horizon = float(contract["modes"][mode]["recommended_word_horizon_s"])
        # Exact source-word starts are block-aligned near the preferred injection
        # time.  Two seconds safely covers the block alignment and final sample.
        required_end.append(float(inject) + horizon + 2.0)
    stop_s = max(required_end)

    prefix_dir = args.diagnostic_dir.resolve() / "_source_prefixes"
    prefix_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    data_dir = args.data_dir.resolve()
    for _, _, name in BASE.RECORDS:
        source = data_dir / name
        if not source.exists():
            raise FileNotFoundError(source)
        manifest.append(write_prefix(source, prefix_dir / name, stop_s))
    (prefix_dir / "prefix_manifest.json").write_text(
        json.dumps(
            {
                "schema": 1,
                "qualification": "EXACT_SOURCE_PREFIX_RUNTIME_OPTIMIZATION_ONLY",
                "modes": modes,
                "required_stop_s": stop_s,
                "records": manifest,
            },
            indent=2,
            sort_keys=True,
        )
    )

    forwarded = replace_option(sys.argv[1:], "--data-dir", str(prefix_dir))
    proc = subprocess.run([sys.executable, str(DRIVER), *forwarded], cwd=REPO)
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
