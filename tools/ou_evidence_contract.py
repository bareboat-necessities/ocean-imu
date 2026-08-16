#!/usr/bin/env python3
"""Normalize OU statistical evidence and bind it to the estimator implementation.

The Monte Carlo studies score continuous metrics and paired effects. Historical
simulator regression thresholds are intentionally excluded from the statistical
bundle so a randomized replay cannot be mistaken for a failed experiment.

The manifest also records a transitive hash closure of repository-local C/C++
implementation dependencies reachable from the simulator targets, plus the
Python analysis/pipeline sources. This makes implementation changes invalidate
committed evidence until the study is regenerated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
GATE_FIELDS = {"quality_gate_pass", "simulator_return_code"}
INCLUDE_RE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.MULTILINE)

STUDIES = {
    "validation": {
        "dir": REPO_ROOT / "reports" / "results" / "ou_validation",
        "raw_csv": "ou_validation_raw.csv",
        "json": "ou_validation.json",
        "manifest": "ou_validation_manifest.json",
        "macro": "ou_validation_macros.tex",
        "mirrored_macro": REPO_ROOT / "doc" / "kalman_ou_iii" / "w3d-ou-validation-macros-generated.tex-part",
        "cpp_roots": (
            REPO_ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim.cpp",
            REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim.cpp",
            REPO_ROOT / "src" / "util" / "W3dSimCommon.cpp",
        ),
        "analysis_files": (
            REPO_ROOT / "tools" / "ou_validation.py",
        ),
    },
    "robustness": {
        "dir": REPO_ROOT / "reports" / "results" / "ou_robustness",
        "raw_csv": "ou_robustness_raw.csv",
        "json": "ou_robustness.json",
        "manifest": "ou_robustness_manifest.json",
        "macro": None,
        "mirrored_macro": None,
        "cpp_roots": (
            REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim.cpp",
            REPO_ROOT / "src" / "util" / "W3dSimCommon.cpp",
        ),
        "analysis_files": (
            REPO_ROOT / "tools" / "ou_validation.py",
            REPO_ROOT / "tools" / "ou_robustness.py",
        ),
    },
}


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {"sha256": sha256_file(path), "bytes": path.stat().st_size}


def repo_name(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def resolve_local_include(source: Path, include: str) -> Path | None:
    candidates = (
        source.parent / include,
        REPO_ROOT / "src" / include,
        REPO_ROOT / include,
    )
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
            resolved.relative_to(REPO_ROOT)
        except (OSError, ValueError):
            continue
        if resolved.is_file():
            return resolved
    return None


def implementation_closure(roots: Iterable[Path]) -> list[Path]:
    pending = [path.resolve() for path in roots]
    seen: set[Path] = set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        if not path.is_file():
            raise FileNotFoundError(path)
        seen.add(path)
        text = path.read_text(encoding="utf-8", errors="ignore")
        for include in INCLUDE_RE.findall(text):
            resolved = resolve_local_include(path, include)
            if resolved is not None and resolved not in seen:
                pending.append(resolved)
    return sorted(seen, key=repo_name)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def strip_gate_csv(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path}")
        fields = [name for name in reader.fieldnames if name not in GATE_FIELDS]
        rows = [
            {name: value for name, value in row.items() if name in fields}
            for row in reader
        ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def strip_gate_json(path: Path) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload.get("raw_runs", []):
        for field in GATE_FIELDS:
            row.pop(field, None)
    protocol = payload.setdefault("protocol", {})
    protocol.pop("quality_gate_window_sec", None)
    protocol["simulator_regression_gates_exported"] = False
    protocol["replay_inclusion_rule"] = (
        "all completed replays with machine-readable metrics; deterministic "
        "simulator regression thresholds are not statistical acceptance criteria"
    )
    write_json(path, payload)


def strip_gate_macros(path: Path | None) -> None:
    if path is None or not path.exists():
        return
    lines = path.read_text(encoding="utf-8").splitlines()
    kept = [
        line for line in lines
        if "OUValidationGatePasses" not in line and "OUValidationGateFailures" not in line
    ]
    path.write_text("\n".join(kept) + "\n", encoding="utf-8")


def stamp_manifest(study: str) -> None:
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    manifest_path = out / str(cfg["manifest"])
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    closure = implementation_closure(cfg["cpp_roots"])
    analysis_files = [Path(path).resolve() for path in cfg["analysis_files"]]
    for path in analysis_files:
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest["implementation_files"] = {
        repo_name(path): file_record(path) for path in closure
    }
    manifest["analysis_pipeline_files"] = {
        repo_name(path): file_record(path) for path in sorted(analysis_files, key=repo_name)
    }
    protocol = manifest.setdefault("protocol", {})
    protocol.pop("quality_gate_window_sec", None)
    protocol["simulator_regression_gates_exported"] = False
    protocol["replay_inclusion_rule"] = (
        "all completed replays with machine-readable metrics; deterministic "
        "simulator regression thresholds are not statistical acceptance criteria"
    )

    result_files: dict[str, dict[str, object]] = {}
    for name in manifest.get("result_files", {}):
        path = out / name
        if path.is_file():
            result_files[name] = file_record(path)
    manifest["result_files"] = result_files
    write_json(manifest_path, manifest)


def normalize(study: str) -> None:
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    strip_gate_csv(out / str(cfg["raw_csv"]))
    strip_gate_json(out / str(cfg["json"]))
    strip_gate_macros(out / str(cfg["macro"]) if cfg["macro"] else None)
    strip_gate_macros(Path(cfg["mirrored_macro"]) if cfg["mirrored_macro"] else None)
    stamp_manifest(study)


def verify_hash_map(records: dict[str, dict[str, object]]) -> list[str]:
    errors: list[str] = []
    for name, expected in records.items():
        path = REPO_ROOT / name
        if not path.is_file():
            errors.append(f"missing provenance file: {name}")
            continue
        current = file_record(path)
        if current != expected:
            errors.append(f"stale provenance: {name}")
    return errors


def check(study: str) -> list[str]:
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    errors: list[str] = []

    raw = out / str(cfg["raw_csv"])
    with raw.open(newline="", encoding="utf-8") as handle:
        fields = csv.DictReader(handle).fieldnames or []
    leaked = GATE_FIELDS.intersection(fields)
    if leaked:
        errors.append(f"{study}: legacy gate columns present: {sorted(leaked)}")

    bundle = json.loads((out / str(cfg["json"])).read_text(encoding="utf-8"))
    if any(GATE_FIELDS.intersection(row) for row in bundle.get("raw_runs", [])):
        errors.append(f"{study}: legacy gate fields present in raw_runs")
    protocol = bundle.get("protocol", {})
    if protocol.get("simulator_regression_gates_exported") is not False:
        errors.append(f"{study}: statistical evidence gate policy is not explicit")

    macro_name = cfg["macro"]
    if macro_name:
        text = (out / str(macro_name)).read_text(encoding="utf-8")
        if "OUValidationGatePasses" in text or "OUValidationGateFailures" in text:
            errors.append("validation: legacy gate-count macros present")

    manifest = json.loads((out / str(cfg["manifest"])).read_text(encoding="utf-8"))
    implementation = manifest.get("implementation_files")
    analysis_files = manifest.get("analysis_pipeline_files")
    if not implementation:
        errors.append(f"{study}: implementation_files missing from manifest")
    else:
        errors.extend(f"{study}: {msg}" for msg in verify_hash_map(implementation))
    if not analysis_files:
        errors.append(f"{study}: analysis_pipeline_files missing from manifest")
    else:
        errors.extend(f"{study}: {msg}" for msg in verify_hash_map(analysis_files))
    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("validation", "robustness", "all"), default="all")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="verify without modifying files")
    mode.add_argument("--auto", action="store_true", help="normalize only freshly regenerated bundles; otherwise verify")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    studies = tuple(STUDIES) if args.study == "all" else (args.study,)
    if args.auto:
        head = git_head()
        fresh = True
        for study in studies:
            cfg = STUDIES[study]
            manifest = json.loads((Path(cfg["dir"]) / str(cfg["manifest"])).read_text(encoding="utf-8"))
            if manifest.get("git_commit") != head:
                fresh = False
                break
        if fresh:
            for study in studies:
                normalize(study)
                errors = check(study)
                if errors:
                    raise RuntimeError("; ".join(errors))
                print(f"Normalized freshly regenerated {study} evidence")
            return 0
        args.check = True

    if args.check:
        errors = [error for study in studies for error in check(study)]
        if errors:
            for error in errors:
                print(error)
            return 1
        print("OU evidence contract is current")
        return 0

    for study in studies:
        normalize(study)
        errors = check(study)
        if errors:
            raise RuntimeError("; ".join(errors))
        print(f"Normalized {study} evidence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
