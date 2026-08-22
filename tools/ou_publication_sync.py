#!/usr/bin/env python3
"""Keep generated OU publication text aligned with committed evidence.

The statistical generator deliberately owns numbers and tables. Editorial
claim scope is owned by the article. This helper performs only mechanical
publication synchronization: it removes one retired interpretation sentence,
refreshes that file's manifest hash, and mirrors the direction tables from the
current validation summary and deterministic OU--III table. It never changes
replay rows or statistics.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


RETIRED_AXIS_CLAIM = (
    " The vertical gain of OU--III is paid for in the horizontal channels."
)
CURRENT_CAPTION_TAIL = r"positive $\Delta$ favors OU--II.}"
PUBLICATION_NAME = "ou_validation_publication.tex"
MANIFEST_NAME = "ou_validation_manifest.json"
SUMMARY_NAME = "ou_validation_summary.csv"

DIRECTION_SCENARIOS = (
    ("stationary_jonswap_H0_270_L14_047_A30_00_P60_00", "JONSWAP", "0.27"),
    ("stationary_jonswap_H1_500_L50_710_A_30_00_P120_00", "JONSWAP", "1.50"),
    ("stationary_jonswap_H4_000_L112_766_A30_00_P30_00", "JONSWAP", "4.00"),
    ("stationary_jonswap_H8_500_L202_839_A_30_00_P72_00", "JONSWAP", "8.50"),
    ("stationary_pmstokes_H0_270_L14_047_A30_00_P60_00", "PM--Stokes", "0.27"),
    ("stationary_pmstokes_H1_500_L50_710_A_30_00_P120_00", "PM--Stokes", "1.50"),
    ("stationary_pmstokes_H4_000_L112_766_A30_00_P30_00", "PM--Stokes", "4.00"),
    ("stationary_pmstokes_H8_500_L202_839_A_30_00_P72_00", "PM--Stokes", "8.50"),
)

DETERMINISTIC_ROW_RE = re.compile(
    r"^\s*(JONSWAP|PM--Stokes)\s*&\s*([0-9.]+)\s*&\s*"
    r"[0-9.]+\s*&\s*[0-9.]+\s*&\s*"
    r"([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*"
    r"(-?[0-9.]+)\s*&\s*([0-9.]+/[0-9.]+/[0-9.]+)\s*\\\\$",
    re.MULTILINE,
)


def _file_record(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def sync_validation_publication(validation_dir: Path) -> bool:
    publication = validation_dir / PUBLICATION_NAME
    manifest_path = validation_dir / MANIFEST_NAME
    if not publication.is_file():
        raise FileNotFoundError(publication)
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)

    text = publication.read_text(encoding="utf-8")
    count = text.count(RETIRED_AXIS_CLAIM)
    if count > 1:
        raise RuntimeError(
            f"retired publication claim occurs {count} times; refusing ambiguous rewrite"
        )
    if count == 1:
        text = text.replace(RETIRED_AXIS_CLAIM, "", 1)
        publication.write_text(text, encoding="utf-8")
    if CURRENT_CAPTION_TAIL not in text:
        raise RuntimeError(
            "OU 3-D publication caption no longer matches the article wording contract"
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    result_files = manifest.get("result_files")
    if not isinstance(result_files, dict) or PUBLICATION_NAME not in result_files:
        raise RuntimeError(
            f"manifest result_files has no {PUBLICATION_NAME} record"
        )
    new_record = _file_record(publication)
    changed = result_files[PUBLICATION_NAME] != new_record or count == 1
    result_files[PUBLICATION_NAME] = new_record
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return changed


def _replace_table_rows(text: str, label: str, rows: list[str]) -> str:
    marker = rf"\label{{{label}}}"
    try:
        label_pos = text.index(marker)
        midrule_end = text.index(r"\midrule", label_pos) + len(r"\midrule")
        bottomrule = text.index(r"\bottomrule", midrule_end)
    except ValueError as exc:
        raise RuntimeError(f"cannot locate table {label}") from exc
    return text[:midrule_end] + "\n" + "\n".join(rows) + "\n    " + text[bottomrule:]


def _direction_summary_rows(validation_dir: Path) -> list[str]:
    summary_path = validation_dir / SUMMARY_NAME
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    wanted = {scenario for scenario, _, _ in DIRECTION_SCENARIOS}
    metrics: dict[tuple[str, str], tuple[float, float]] = {}
    with summary_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            scenario = row.get("scenario", "")
            metric = row.get("metric", "")
            if (
                scenario in wanted
                and row.get("family") == "OU_III"
                and row.get("mode") == "Adaptive"
                and metric in {"dir_axis_rmse_deg", "dir_travel_rmse_deg"}
            ):
                metrics[(scenario, metric)] = (
                    float(row["mean"]),
                    float(row["std"]),
                )

    expected = {
        (scenario, metric)
        for scenario, _, _ in DIRECTION_SCENARIOS
        for metric in ("dir_axis_rmse_deg", "dir_travel_rmse_deg")
    }
    missing = sorted(expected - set(metrics))
    if missing:
        raise RuntimeError(f"validation summary is missing direction rows: {missing}")

    rows: list[str] = []
    for index, (scenario, case, hs) in enumerate(DIRECTION_SCENARIOS):
        axis_mean, axis_std = metrics[(scenario, "dir_axis_rmse_deg")]
        travel_mean, travel_std = metrics[(scenario, "dir_travel_rmse_deg")]
        rows.append(
            f"    {case:<11} & {hs} & "
            f"${axis_mean:.2f}\\pm{axis_std:.2f}$ & "
            f"${travel_mean:.2f}\\pm{travel_std:.2f}$ \\\\"
        )
        if index == 3:
            rows.append(r"    \addlinespace")
    return rows


def _deterministic_direction_rows(deterministic_results: Path) -> list[str]:
    if not deterministic_results.is_file():
        raise FileNotFoundError(deterministic_results)
    matches = DETERMINISTIC_ROW_RE.findall(
        deterministic_results.read_text(encoding="utf-8")
    )
    if len(matches) != 8:
        raise RuntimeError(
            f"expected 8 deterministic OU--III rows, found {len(matches)}"
        )

    rows: list[str] = []
    for index, (case, hs, roll, pitch, yaw, theta, tau) in enumerate(matches):
        rows.append(
            f"    {case:<11} & {hs} & {roll} & {pitch} & {yaw} & "
            f"{theta:>5} & {tau} \\\\"
        )
        if index == 3:
            rows.append(r"    \addlinespace")
    return rows


def sync_direction_tables(
    validation_dir: Path,
    direction_results: Path,
    deterministic_results: Path,
) -> bool:
    if not direction_results.is_file():
        raise FileNotFoundError(direction_results)

    original = direction_results.read_text(encoding="utf-8")
    text = _replace_table_rows(
        original,
        "tab:direction-ou3-rms",
        _direction_summary_rows(validation_dir),
    )
    text = _replace_table_rows(
        text,
        "tab:direction-ou3-integration",
        _deterministic_direction_rows(deterministic_results),
    )
    if text == original:
        return False
    direction_results.write_text(text, encoding="utf-8")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, required=True)
    parser.add_argument("--direction-results", type=Path)
    parser.add_argument("--deterministic-results", type=Path)
    args = parser.parse_args()

    if bool(args.direction_results) != bool(args.deterministic_results):
        parser.error(
            "--direction-results and --deterministic-results must be supplied together"
        )

    publication_changed = sync_validation_publication(args.validation_dir)
    direction_changed = False
    if args.direction_results is not None:
        direction_changed = sync_direction_tables(
            args.validation_dir,
            args.direction_results,
            args.deterministic_results,
        )

    if publication_changed or direction_changed:
        print("Aligned OU publication wording and direction tables with evidence.")
    else:
        print("OU publication wording and direction tables already aligned.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
