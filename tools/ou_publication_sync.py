#!/usr/bin/env python3
"""Synchronize OU publication inputs while enforcing current claim scope."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

RETIRED_AXIS_CLAIM = " The vertical gain of OU--III is paid for in the horizontal channels."
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
VALIDATION_ARTICLE_SVGS = (
    "ou_validation_vertical.svg",
    "ou_validation_displacement.svg",
    "ou_validation_attitude.svg",
)
ROBUSTNESS_PUBLICATION_NAME = "ou_robustness_publication.tex"
ROBUSTNESS_MACROS_NAME = "ou_robustness_macros.tex"
ROBUSTNESS_STRESS_SVG = "ou_robustness_stress.svg"
ROBUSTNESS_DOC_RESULTS = "w3d-ou-robustness-results-generated.tex-part"
ROBUSTNESS_DOC_MACROS = "w3d-ou-robustness-macros-generated.tex-part"
ROBUSTNESS_RETIRED_SVG = "ou_robustness_sensitivity.svg"
ROBUSTNESS_STRESS_MACROS = (
    "OURobustnessPairs", "OURobustnessLowReferenceMean", "OURobustnessLowStressMean",
    "OURobustnessLowStressAbsolute", "OURobustnessLowDifference",
    "OURobustnessLowDifferenceLow", "OURobustnessLowDifferenceHigh",
    "OURobustnessControlledAdaptiveMean", "OURobustnessRapidAdaptiveMean",
    "OURobustnessRapidDifference", "OURobustnessRapidDifferenceLow",
    "OURobustnessRapidDifferenceHigh", "OURobustnessRapidAdaptationDifference",
    "OURobustnessRapidAdaptationLow", "OURobustnessRapidAdaptationHigh",
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
    r"^\s*(JONSWAP|PM--Stokes)\s*&\s*([0-9.]+)\s*&\s*[0-9.]+\s*&\s*[0-9.]+\s*&\s*"
    r"([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*(-?[0-9.]+)\s*&\s*"
    r"([0-9.]+/[0-9.]+/[0-9.]+)\s*\\\\$",
    re.MULTILINE,
)


def _file_record(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def _replace_caption_before_label(text: str, label: str, caption: str) -> tuple[str, bool]:
    marker = rf"\label{{{label}}}"
    pos = text.index(marker)
    cap = text.rfind(r"\caption{", 0, pos)
    line = text.rfind("\n", 0, pos) + 1
    current = text[cap:line].rstrip()
    if current == caption:
        return text, False
    return text[:cap] + caption + "\n" + text[line:], True


def _strip_table(text: str, label: str) -> str:
    marker = rf"\label{{{label}}}"
    pos = text.find(marker)
    if pos < 0:
        return text
    start = text.rfind(r"\begin{table*}", 0, pos)
    token = r"\end{table*}"
    end = text.find(token, pos)
    if start < 0 or end < 0:
        raise RuntimeError(f"cannot delimit table {label}")
    return text[:start].rstrip() + "\n\n" + text[end + len(token):].lstrip()


def _replace_table_caption(text: str, label: str, caption: str) -> str:
    marker = rf"\label{{{label}}}"
    pos = text.find(marker)
    if pos < 0:
        return text
    start = text.rfind(r"\begin{table*}", 0, pos)
    cap = text.find(r"\caption{", start, pos)
    if start < 0 or cap < 0:
        raise RuntimeError(f"caption not found for table {label}")
    body_start = cap + len(r"\caption{")
    depth = 1
    i = body_start
    while i < pos and depth:
        if text[i] == "{" and (i == 0 or text[i - 1] != "\\"):
            depth += 1
        elif text[i] == "}" and (i == 0 or text[i - 1] != "\\"):
            depth -= 1
        i += 1
    if depth != 0:
        raise RuntimeError(f"unbalanced caption for table {label}")
    return text[:cap] + r"\caption{" + caption + "}" + text[i:]


def _compact_axes_table(text: str) -> str:
    label = r"\label{tab:ou_mc_axes}"
    pos = text.find(label)
    if pos < 0:
        raise RuntimeError("3-D axis table label not found")
    start = text.rfind(r"\begin{table*}", 0, pos)
    end_token = r"\end{table*}"
    end = text.find(end_token, pos)
    if start < 0 or end < 0:
        raise RuntimeError("could not delimit 3-D axis table")
    end += len(end_token)
    block = text[start:end]

    spacing = r"\setlength{\tabcolsep}{3.0pt}"
    if block.count(spacing) == 1:
        block = block.replace(spacing, r"\setlength{\tabcolsep}{1.6pt}", 1)
    elif r"\setlength{\tabcolsep}{1.6pt}" not in block:
        raise RuntimeError("3-D axis table spacing anchor moved")

    block, scenario_count = re.subn(r"\$H_s=([0-9.]+)\$ m", r"\1 m", block)
    if scenario_count not in (0, 4):
        raise RuntimeError(
            f"expected zero or four stationary H_s labels in 3-D axis table, found {scenario_count}"
        )
    block = re.sub(
        r"([0-9]+\.[0-9]+) \$\\pm\$ ([0-9]+\.[0-9]+)",
        r"$\1{\\pm}\2$",
        block,
    )
    return text[:start] + block + text[end:]


def curate_validation_for_article(text: str) -> str:
    """Create the compact article view and exclude the retired one-way transition."""
    for label in (
        "tab:ou_mc_adaptation",
        "tab:ou_mc_covsync",
        "tab:ou_mc_direction",
        "tab:ou_transition_segments",
    ):
        text = _strip_table(text, label)

    text = re.sub(r"(?m)^\s*Transition\s*&.*?\\\\\s*$\n?", "", text)
    text = _compact_axes_table(text)
    for label, caption in (
        ("tab:ou_mc_family", r"Ten-seed paired OU-family comparison over the final \SI{900}{s}."),
        ("tab:ou_mc_axes", r"Ten-seed adaptive displacement RMS by axis over the final \SI{900}{s}."),
        ("tab:ou_mc_channels", r"OU--III adaptation-channel ablation over the final \SI{900}{s}."),
        ("tab:ou_mc_pmstokes", r"Ten-seed paired OU-family comparison on PM--Stokes seas."),
    ):
        text = _replace_table_caption(text, label, caption)

    if "tab:ou_transition_segments" in text or re.search(r"(?m)^\s*Transition\s*&", text):
        raise RuntimeError("one-way transition survived publication curation")
    return text


def sync_validation_publication(validation_dir: Path) -> bool:
    publication = validation_dir / PUBLICATION_NAME
    manifest_path = validation_dir / MANIFEST_NAME
    text = publication.read_text(encoding="utf-8")
    changed = False
    count = text.count(RETIRED_AXIS_CLAIM)
    if count > 1:
        raise RuntimeError("retired axis claim occurs more than once")
    if count:
        text = text.replace(RETIRED_AXIS_CLAIM, "", 1)
        changed = True
    if CURRENT_CAPTION_TAIL not in text:
        raise RuntimeError("validation caption contract moved")
    text, caption_changed = _replace_caption_before_label(
        text, "tab:ou_mc_channels", CURRENT_CHANNEL_CAPTION
    )
    changed |= caption_changed
    if changed:
        publication.write_text(text, encoding="utf-8")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    files = manifest["result_files"]
    record = _file_record(publication)
    if files[PUBLICATION_NAME] != record:
        files[PUBLICATION_NAME] = record
        changed = True
    if changed:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return changed


def sync_validation_doc_copy(validation_dir: Path, doc_dir: Path) -> bool:
    source = (validation_dir / PUBLICATION_NAME).read_text(encoding="utf-8")
    changed = False
    generated = doc_dir / "w3d-ou-validation-results-generated.tex-part"
    if not generated.is_file() or generated.read_text(encoding="utf-8") != source:
        generated.write_text(source, encoding="utf-8")
        changed = True
    publication = doc_dir / "w3d-ou-validation-results-publication.tex-part"
    desired = curate_validation_for_article(source)
    if not publication.is_file() or publication.read_text(encoding="utf-8") != desired:
        publication.write_text(desired, encoding="utf-8")
        changed = True
    for name in VALIDATION_ARTICLE_SVGS:
        source_svg = validation_dir / name
        target_svg = doc_dir / name
        svg = source_svg.read_bytes()
        if not target_svg.is_file() or target_svg.read_bytes() != svg:
            target_svg.write_bytes(svg)
            changed = True
    retired = doc_dir / "ou_validation_transition.svg"
    if retired.exists():
        retired.unlink()
        changed = True
    return changed


def _table_block_by_label(text: str, label: str) -> str:
    marker = rf"\label{{{label}}}"
    for block in re.findall(r"\\begin\{table\*\}.*?\\end\{table\*\}", text, flags=re.S):
        if marker in block:
            return block
    raise RuntimeError(f"cannot locate table {label}")


def _robustness_macro_definitions(text: str) -> dict[str, str]:
    return dict(re.findall(r"\\providecommand\{\\(OURobustness[A-Za-z]+)\}\{([^}]*)\}", text))


def sync_robustness_doc_copies(robustness_dir: Path, doc_dir: Path) -> bool:
    publication = robustness_dir / ROBUSTNESS_PUBLICATION_NAME
    macros_path = robustness_dir / ROBUSTNESS_MACROS_NAME
    stress_svg = robustness_dir / ROBUSTNESS_STRESS_SVG
    changed = False
    table = _table_block_by_label(
        publication.read_text(encoding="utf-8"), "tab:ou_robustness_stress"
    )
    desired = "% Publication excerpt of the committed OU--III degradation cases.\n\n" + table + "\n"
    target = doc_dir / ROBUSTNESS_DOC_RESULTS
    if not target.is_file() or target.read_text(encoding="utf-8") != desired:
        target.write_text(desired, encoding="utf-8")
        changed = True
    macros = _robustness_macro_definitions(macros_path.read_text(encoding="utf-8"))
    missing = [name for name in ROBUSTNESS_STRESS_MACROS if name not in macros]
    if missing:
        raise RuntimeError(f"robustness archive is missing macros: {missing}")
    desired_macros = (
        "% Publication subset of macros generated from the committed OU--III robustness bundle.\n"
        + "".join(
            f"\\providecommand{{\\{name}}}{{{macros[name]}}}\n"
            for name in ROBUSTNESS_STRESS_MACROS
        )
    )
    macros_target = doc_dir / ROBUSTNESS_DOC_MACROS
    if not macros_target.is_file() or macros_target.read_text(encoding="utf-8") != desired_macros:
        macros_target.write_text(desired_macros, encoding="utf-8")
        changed = True
    stress_target = doc_dir / ROBUSTNESS_STRESS_SVG
    svg = stress_svg.read_bytes()
    if not stress_target.is_file() or stress_target.read_bytes() != svg:
        stress_target.write_bytes(svg)
        changed = True
    retired = doc_dir / ROBUSTNESS_RETIRED_SVG
    if retired.exists():
        retired.unlink()
        changed = True
    return changed


def _replace_table_rows(text: str, label: str, rows: list[str]) -> str:
    pos = text.index(rf"\label{{{label}}}")
    mid = text.index(r"\midrule", pos) + len(r"\midrule")
    bottom = text.index(r"\bottomrule", mid)
    return text[:mid] + "\n" + "\n".join(rows) + "\n    " + text[bottom:]


def _direction_summary_rows(validation_dir: Path) -> list[str]:
    wanted = {scenario for scenario, _, _ in DIRECTION_SCENARIOS}
    metrics = {}
    with (validation_dir / SUMMARY_NAME).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = (row.get("scenario", ""), row.get("metric", ""))
            if (
                row.get("scenario") in wanted
                and row.get("family") == "OU_III"
                and row.get("mode") == "Adaptive"
                and row.get("metric") in {"dir_axis_rmse_deg", "dir_travel_rmse_deg"}
            ):
                metrics[key] = (float(row["mean"]), float(row["std"]))
    rows = []
    for index, (scenario, case, hs) in enumerate(DIRECTION_SCENARIOS):
        axis, axis_std = metrics[(scenario, "dir_axis_rmse_deg")]
        travel, travel_std = metrics[(scenario, "dir_travel_rmse_deg")]
        rows.append(
            f"    {case:<11} & {hs} & ${axis:.2f}\\pm{axis_std:.2f}$ & "
            f"${travel:.2f}\\pm{travel_std:.2f}$ \\\\"
        )
        if index == 3:
            rows.append(r"    \addlinespace")
    return rows


def _deterministic_direction_rows(path: Path) -> list[str]:
    matches = DETERMINISTIC_ROW_RE.findall(path.read_text(encoding="utf-8"))
    if len(matches) != 8:
        raise RuntimeError(f"expected 8 deterministic rows, found {len(matches)}")
    rows = []
    for index, (case, hs, roll, pitch, yaw, theta, tau) in enumerate(matches):
        rows.append(
            f"    {case:<11} & {hs} & {roll} & {pitch} & {yaw} & {theta:>5} & {tau} \\\\"
        )
        if index == 3:
            rows.append(r"    \addlinespace")
    return rows


def sync_direction_tables(validation_dir: Path, direction_results: Path, deterministic_results: Path) -> bool:
    original = direction_results.read_text(encoding="utf-8")
    text = _replace_table_rows(
        original, "tab:direction-ou3-rms", _direction_summary_rows(validation_dir)
    )
    text = _replace_table_rows(
        text, "tab:direction-ou3-integration", _deterministic_direction_rows(deterministic_results)
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
        parser.error("direction arguments must be supplied together")
    root = Path(__file__).resolve().parents[1]
    doc = root / "doc" / "kalman_ou_iii"
    changed = sync_validation_publication(args.validation_dir) | sync_validation_doc_copy(args.validation_dir, doc)
    robustness = args.validation_dir.parent / "ou_robustness"
    if robustness.is_dir():
        changed |= sync_robustness_doc_copies(robustness, doc)
    if args.direction_results is not None:
        changed |= sync_direction_tables(
            args.validation_dir, args.direction_results, args.deterministic_results
        )
    print(
        "Aligned OU publication tree with current evidence and claim scope."
        if changed
        else "OU publication tree already aligned with current evidence and claim scope."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
