#!/usr/bin/env python3
"""Keep OU publication text aligned with evidence without republishing retired studies.

The statistical generators own replay rows, statistics, and the complete evidence
archives. The article owns claim scope. This helper aligns editable validation
wording, mirrors the current validation publication into the manuscript, and
curates only law-independent degradation evidence from the robustness archive
into the OU--III publication tree. Historical robustness sensitivity evidence
remains under reports/results for provenance.
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
CURRENT_CHANNEL_CAPTION = (
    r"\caption{OU--III adaptation-channel ablation for vertical-displacement RMS "
    r"error over the final \SI{900}{s}, in percent of $H_s$ (mean $\pm$ sample "
    r"standard deviation, $n=10$ paired seed triplets). The four columns form a "
    r"$2\times2$ factorial in the applied parameter channels. \emph{$r_S$ only} "
    r"freezes $\tau$ and $\sigma_{aw}$ while the deployed SpectralMSE regularizer "
    r"channel continues to adapt; \emph{OU only} adapts $\tau$ and $\sigma_{aw}$ "
    r"while holding $r_S$ at FixedNominal. This isolates the two applied adaptation "
    r"channels without introducing or comparing an alternative regularizer law.}"
)
PUBLICATION_NAME = "ou_validation_publication.tex"
MANIFEST_NAME = "ou_validation_manifest.json"
SUMMARY_NAME = "ou_validation_summary.csv"

ROBUSTNESS_PUBLICATION_NAME = "ou_robustness_publication.tex"
ROBUSTNESS_MACROS_NAME = "ou_robustness_macros.tex"
ROBUSTNESS_STRESS_SVG = "ou_robustness_stress.svg"
ROBUSTNESS_DOC_RESULTS = "w3d-ou-robustness-results-generated.tex-part"
ROBUSTNESS_DOC_MACROS = "w3d-ou-robustness-macros-generated.tex-part"
ROBUSTNESS_RETIRED_SVG = "ou_robustness_sensitivity.svg"
ROBUSTNESS_STRESS_MACROS = (
    "OURobustnessPairs",
    "OURobustnessLowReferenceMean",
    "OURobustnessLowStressMean",
    "OURobustnessLowStressAbsolute",
    "OURobustnessLowDifference",
    "OURobustnessLowDifferenceLow",
    "OURobustnessLowDifferenceHigh",
    "OURobustnessControlledAdaptiveMean",
    "OURobustnessRapidAdaptiveMean",
    "OURobustnessRapidDifference",
    "OURobustnessRapidDifferenceLow",
    "OURobustnessRapidDifferenceHigh",
    "OURobustnessRapidAdaptationDifference",
    "OURobustnessRapidAdaptationLow",
    "OURobustnessRapidAdaptationHigh",
)

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


def _replace_caption_before_label(
    text: str, label: str, caption: str
) -> tuple[str, bool]:
    marker = rf"\label{{{label}}}"
    try:
        label_pos = text.index(marker)
    except ValueError as exc:
        raise RuntimeError(f"cannot locate table label {label}") from exc
    caption_pos = text.rfind(r"\caption{", 0, label_pos)
    if caption_pos < 0:
        raise RuntimeError(f"cannot locate caption before table {label}")
    label_line_start = text.rfind("\n", 0, label_pos) + 1
    current = text[caption_pos:label_line_start].rstrip()
    if current == caption:
        return text, False
    return text[:caption_pos] + caption + "\n" + text[label_line_start:], True


def sync_validation_publication(validation_dir: Path) -> bool:
    publication = validation_dir / PUBLICATION_NAME
    manifest_path = validation_dir / MANIFEST_NAME
    if not publication.is_file():
        raise FileNotFoundError(publication)
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)

    text = publication.read_text(encoding="utf-8")
    changed = False

    count = text.count(RETIRED_AXIS_CLAIM)
    if count > 1:
        raise RuntimeError(
            f"retired publication claim occurs {count} times; refusing ambiguous rewrite"
        )
    if count == 1:
        text = text.replace(RETIRED_AXIS_CLAIM, "", 1)
        changed = True

    if CURRENT_CAPTION_TAIL not in text:
        raise RuntimeError(
            "OU 3-D publication caption no longer matches the article wording contract"
        )

    text, channel_changed = _replace_caption_before_label(
        text, "tab:ou_mc_channels", CURRENT_CHANNEL_CAPTION
    )
    changed = changed or channel_changed
    if changed:
        publication.write_text(text, encoding="utf-8")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    result_files = manifest.get("result_files")
    if not isinstance(result_files, dict) or PUBLICATION_NAME not in result_files:
        raise RuntimeError(f"manifest result_files has no {PUBLICATION_NAME} record")

    new_record = _file_record(publication)
    if result_files[PUBLICATION_NAME] != new_record:
        result_files[PUBLICATION_NAME] = new_record
        changed = True
    if changed:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return changed


def sync_validation_doc_copy(validation_dir: Path, doc_dir: Path) -> bool:
    source = validation_dir / PUBLICATION_NAME
    target = doc_dir / "w3d-ou-validation-results-generated.tex-part"
    desired = source.read_bytes()
    if target.is_file() and target.read_bytes() == desired:
        return False
    target.write_bytes(desired)
    return True


def _table_block_by_label(text: str, label: str) -> str:
    marker = rf"\label{{{label}}}"
    for block in re.findall(
        r"\\begin\{table\*\}.*?\\end\{table\*\}", text, flags=re.S
    ):
        if marker in block:
            return block
    raise RuntimeError(f"cannot locate table {label}")


def _robustness_macro_definitions(text: str) -> dict[str, str]:
    return dict(
        re.findall(
            r"\\providecommand\{\\(OURobustness[A-Za-z]+)\}\{([^}]*)\}", text
        )
    )


def sync_robustness_doc_copies(robustness_dir: Path, doc_dir: Path) -> bool:
    publication = robustness_dir / ROBUSTNESS_PUBLICATION_NAME
    macros_path = robustness_dir / ROBUSTNESS_MACROS_NAME
    stress_svg = robustness_dir / ROBUSTNESS_STRESS_SVG
    for path in (publication, macros_path, stress_svg):
        if not path.is_file():
            raise FileNotFoundError(path)

    changed = False

    stress_table = _table_block_by_label(
        publication.read_text(encoding="utf-8"), "tab:ou_robustness_stress"
    )
    desired_results = (
        "% Publication excerpt of the committed OU--III degradation cases.\n\n"
        + stress_table
        + "\n"
    )
    results_target = doc_dir / ROBUSTNESS_DOC_RESULTS
    if (
        not results_target.is_file()
        or results_target.read_text(encoding="utf-8") != desired_results
    ):
        results_target.write_text(desired_results, encoding="utf-8")
        changed = True

    archived_macros = _robustness_macro_definitions(
        macros_path.read_text(encoding="utf-8")
    )
    missing = [
        name for name in ROBUSTNESS_STRESS_MACROS if name not in archived_macros
    ]
    if missing:
        raise RuntimeError(f"robustness archive is missing publication macros: {missing}")

    desired_macros = (
        "% Publication subset of macros generated from the committed OU--III "
        "robustness bundle.\n"
        + "".join(
            f"\\providecommand{{\\{name}}}{{{archived_macros[name]}}}\n"
            for name in ROBUSTNESS_STRESS_MACROS
        )
    )
    macros_target = doc_dir / ROBUSTNESS_DOC_MACROS
    if (
        not macros_target.is_file()
        or macros_target.read_text(encoding="utf-8") != desired_macros
    ):
        macros_target.write_text(desired_macros, encoding="utf-8")
        changed = True

    stress_target = doc_dir / ROBUSTNESS_STRESS_SVG
    desired_svg = stress_svg.read_bytes()
    if not stress_target.is_file() or stress_target.read_bytes() != desired_svg:
        stress_target.write_bytes(desired_svg)
        changed = True

    retired_svg = doc_dir / ROBUSTNESS_RETIRED_SVG
    if retired_svg.exists():
        retired_svg.unlink()
        changed = True
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

    repo_root = Path(__file__).resolve().parents[1]
    doc_dir = repo_root / "doc" / "kalman_ou_iii"

    publication_changed = sync_validation_publication(args.validation_dir)
    validation_doc_changed = sync_validation_doc_copy(args.validation_dir, doc_dir)

    robustness_changed = False
    robustness_dir = args.validation_dir.parent / "ou_robustness"
    if robustness_dir.is_dir():
        robustness_changed = sync_robustness_doc_copies(robustness_dir, doc_dir)

    direction_changed = False
    if args.direction_results is not None:
        direction_changed = sync_direction_tables(
            args.validation_dir,
            args.direction_results,
            args.deterministic_results,
        )

    if (
        publication_changed
        or validation_doc_changed
        or robustness_changed
        or direction_changed
    ):
        print("Aligned OU publication tree with current evidence and claim scope.")
    else:
        print("OU publication tree already aligned with current evidence and claim scope.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
