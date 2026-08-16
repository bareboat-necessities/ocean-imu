#!/usr/bin/env python3
"""Shared provenance primitives for OU validation and robustness evidence.

The scientific replay and a later statistical restatement are different events.
`replay_provenance` is immutable once simulator rows exist; only a full simulator
regeneration may create it.  A restatement may update derived statistics and its
own analysis provenance, but it must first prove that the current replay-producing
implementation still matches the replay dependency closure.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import io
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 2
GATE_FIELDS = ("quality_gate_pass", "simulator_return_code")
INCLUDE_RE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.MULTILINE)

STUDIES: dict[str, dict[str, Any]] = {
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
        "build_files": (
            REPO_ROOT / "tests" / "kalman_ou_ii" / "Makefile",
            REPO_ROOT / "tests" / "kalman_ou_iii" / "Makefile",
        ),
        "analysis_files": (REPO_ROOT / "tools" / "ou_validation.py",),
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
        "build_files": (REPO_ROOT / "tests" / "kalman_ou_iii" / "Makefile",),
        "analysis_files": (
            REPO_ROOT / "tools" / "ou_validation.py",
            REPO_ROOT / "tools" / "ou_robustness.py",
        ),
    },
}


class ProvenanceError(RuntimeError):
    pass


@dataclass(frozen=True)
class RestatementContext:
    study: str
    source_bundle: Path
    source_manifest: Path
    source_raw_csv: Path
    replay_provenance: dict[str, Any]
    source_bundle_sha256: str
    source_manifest_sha256: str


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def normalize_record(record: Mapping[str, Any]) -> dict[str, object]:
    size = record.get("size_bytes", record.get("bytes"))
    return {"sha256": str(record.get("sha256", "")), "size_bytes": int(size or 0)}


def legacy_record(record: Mapping[str, Any]) -> dict[str, object]:
    normalized = normalize_record(record)
    return {"sha256": normalized["sha256"], "bytes": normalized["size_bytes"]}


def repo_name(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def resolve_local_include(source: Path, include: str) -> Path | None:
    for candidate in (source.parent / include, REPO_ROOT / "src" / include, REPO_ROOT / include):
        try:
            resolved = candidate.resolve()
            resolved.relative_to(REPO_ROOT)
        except (OSError, ValueError):
            continue
        if resolved.is_file():
            return resolved
    return None


def implementation_closure(roots: Iterable[Path]) -> list[Path]:
    """Repository-local dependency closure.

    Simulator Makefiles already compile with -MMD -MP.  The committed evidence
    contract intentionally keeps this source crawler because combine/archive
    validation may run without compiler-generated .d files.  `build_files` are
    included separately so compiler flags cannot change unnoticed.  The API is
    isolated here so .d-file ingestion can replace the crawler later.
    """
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


def implementation_paths(study: str) -> list[Path]:
    cfg = STUDIES[study]
    paths = set(implementation_closure(cfg["cpp_roots"]))
    paths.update(Path(path).resolve() for path in cfg["build_files"])
    return sorted(paths, key=repo_name)


def implementation_records(study: str) -> dict[str, dict[str, object]]:
    return {repo_name(path): file_record(path) for path in implementation_paths(study)}


def analysis_records(study: str) -> dict[str, dict[str, object]]:
    paths = [Path(path).resolve() for path in STUDIES[study]["analysis_files"]]
    return {repo_name(path): file_record(path) for path in sorted(paths, key=repo_name)}


def git_available() -> bool:
    try:
        completed = subprocess.run(
            ("git", "rev-parse", "--is-inside-work-tree"), cwd=REPO_ROOT,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=False,
        )
    except OSError:
        return False
    return completed.returncode == 0 and completed.stdout.strip() == "true"


def git_output(*args: str) -> str | None:
    if not git_available():
        return None
    completed = subprocess.run(
        ("git", *args), cwd=REPO_ROOT, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def git_bytes(commit: str, path: str) -> bytes:
    completed = subprocess.run(
        ("git", "show", f"{commit}:{path}"), cwd=REPO_ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if completed.returncode != 0:
        raise ProvenanceError(
            f"cannot read {path} at replay commit {commit}: "
            + completed.stderr.decode("utf-8", errors="replace").strip()
        )
    return completed.stdout


def _compiler_environment() -> dict[str, Any]:
    compiler = os.environ.get("CXX", "g++")
    try:
        banner = subprocess.run(
            (compiler, "--version"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, check=False,
        ).stdout.strip()
        version = subprocess.run(
            (compiler, "-dumpfullversion", "-dumpversion"), stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True, check=False,
        ).stdout.strip()
    except OSError:
        banner, version = "unavailable", "unavailable"
    return {
        "command": compiler,
        "version": version or "unknown",
        "banner": banner.splitlines()[0] if banner else "unknown",
        "CXXFLAGS_env": os.environ.get("CXXFLAGS", ""),
        "CPPFLAGS_env": os.environ.get("CPPFLAGS", ""),
    }


def _eigen_environment() -> dict[str, Any]:
    candidates = (
        REPO_ROOT / "third_party" / "eigen" / "Eigen" / "src" / "Core" / "util" / "Macros.h",
        Path("/usr/include/eigen3/Eigen/src/Core/util/Macros.h"),
    )
    for path in candidates:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        values: dict[str, str] = {}
        for macro in ("EIGEN_WORLD_VERSION", "EIGEN_MAJOR_VERSION", "EIGEN_MINOR_VERSION"):
            match = re.search(rf"^\s*#\s*define\s+{macro}\s+(\d+)", text, re.MULTILINE)
            if match:
                values[macro] = match.group(1)
        version = ".".join(
            values.get(name, "?")
            for name in ("EIGEN_WORLD_VERSION", "EIGEN_MAJOR_VERSION", "EIGEN_MINOR_VERSION")
        )
        return {"version": version, "identity_file_sha256": sha256_file(path), "path": str(path)}
    return {"version": "unavailable", "path": "unavailable"}


def environment_metadata() -> dict[str, Any]:
    try:
        import numpy as np  # type: ignore
        numpy_version = np.__version__
    except Exception:
        numpy_version = "unavailable"
    return {
        "python": sys.version,
        "numpy": numpy_version,
        "platform": platform.platform(),
        "compiler": _compiler_environment(),
        "eigen": _eigen_environment(),
    }


def workflow_metadata() -> dict[str, str] | None:
    run_id = os.environ.get("GITHUB_RUN_ID")
    if not run_id:
        return None
    return {
        "run_id": run_id,
        "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT", ""),
        "workflow": os.environ.get("GITHUB_WORKFLOW", ""),
        "job": os.environ.get("GITHUB_JOB", ""),
    }


def _normalized_csv_bytes(data: bytes) -> bytes:
    text = data.decode("utf-8")
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames is None:
        raise ProvenanceError("raw replay CSV has no header")
    fields = [name for name in reader.fieldnames if name not in GATE_FIELDS]
    rows = [{name: row.get(name, "") for name in fields} for row in reader]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def normalized_raw_rows_equal_historical(current_raw: Path, historical_raw: bytes) -> bool:
    current = current_raw.read_bytes()
    return current == _normalized_csv_bytes(historical_raw)


def strip_gate_csv(path: Path) -> None:
    normalized = _normalized_csv_bytes(path.read_bytes())
    path.write_bytes(normalized)


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
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def strip_gate_macros(path: Path | None) -> None:
    if path is None or not path.exists():
        return
    kept = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if "OUValidationGatePasses" not in line and "OUValidationGateFailures" not in line
    ]
    path.write_text("\n".join(kept) + "\n", encoding="utf-8")


def normalize_statistical_schema(study: str) -> None:
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    strip_gate_csv(out / cfg["raw_csv"])
    strip_gate_json(out / cfg["json"])
    strip_gate_macros(out / cfg["macro"] if cfg["macro"] else None)
    strip_gate_macros(Path(cfg["mirrored_macro"]) if cfg["mirrored_macro"] else None)


def _input_records_from_legacy(manifest: Mapping[str, Any]) -> dict[str, dict[str, object]]:
    return {name: normalize_record(record) for name, record in manifest.get("source_files", {}).items()}


def _legacy_map(records: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, object]]:
    return {name: legacy_record(record) for name, record in records.items()}


def _replay_aliases(manifest: dict[str, Any]) -> None:
    replay = manifest["replay_provenance"]
    manifest["git_commit"] = replay.get("git_commit")
    manifest["implementation_files"] = _legacy_map(replay.get("implementation_files", {}))
    manifest["source_files"] = _legacy_map(replay.get("input_files", {}))
    restatement = manifest.get("restatement")
    if restatement and restatement.get("analysis_pipeline_files"):
        manifest["analysis_pipeline_files"] = _legacy_map(restatement["analysis_pipeline_files"])
    manifest.pop("sources_moved_since_rows", None)


def create_replay_provenance_from_fresh_generation(study: str, manifest: Mapping[str, Any]) -> dict[str, Any]:
    if not git_available():
        raise ProvenanceError("full replay provenance creation requires a Git checkout")
    head = git_output("rev-parse", "HEAD")
    recorded = manifest.get("git_commit")
    command = [str(item) for item in manifest.get("command", [])]
    restated = manifest.get("restated_from") or manifest.get("protocol", {}).get("restated_from")
    if not head or recorded != head or restated or "--restat-from" in command:
        raise ProvenanceError(
            "refusing to create replay provenance: the bundle is not a freshly "
            "generated full replay at the current HEAD"
        )
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    raw = out / cfg["raw_csv"]
    return {
        "git_commit": head,
        "implementation_files": implementation_records(study),
        "input_files": _input_records_from_legacy(manifest),
        "raw_rows_sha256": sha256_file(raw),
        "raw_rows_size_bytes": raw.stat().st_size,
        "raw_rows_schema": "normalized statistical replay rows; deterministic simulator gate columns excluded",
        "workflow_run": workflow_metadata(),
        "generated_at": utc_now(),
        "environment": environment_metadata(),
    }


def _compare_record_maps(expected: Mapping[str, Mapping[str, Any]], actual: Mapping[str, Mapping[str, Any]]) -> list[str]:
    errors: list[str] = []
    expected_norm = {name: normalize_record(record) for name, record in expected.items()}
    actual_norm = {name: normalize_record(record) for name, record in actual.items()}
    for name in sorted(set(expected_norm) | set(actual_norm)):
        if name not in expected_norm:
            errors.append(f"new replay dependency not present in replay provenance: {name}")
        elif name not in actual_norm:
            errors.append(f"recorded replay dependency missing from source tree: {name}")
        elif expected_norm[name] != actual_norm[name]:
            errors.append(f"replay dependency differs from replay provenance: {name}")
    return errors


def replay_errors(study: str, manifest: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    replay = manifest.get("replay_provenance")
    if not isinstance(replay, Mapping):
        return [f"{study}: immutable replay_provenance is missing"]
    raw = Path(STUDIES[study]["dir"]) / STUDIES[study]["raw_csv"]
    expected_raw = str(replay.get("raw_rows_sha256", ""))
    if not raw.is_file():
        errors.append(f"{study}: raw replay rows are missing: {raw}")
    elif sha256_file(raw) != expected_raw:
        errors.append(f"{study}: raw replay rows differ from immutable replay provenance")
    expected_impl = replay.get("implementation_files", {})
    if not expected_impl:
        errors.append(f"{study}: replay implementation dependency hashes are missing")
    else:
        errors.extend(f"{study}: {msg}" for msg in _compare_record_maps(expected_impl, implementation_records(study)))
    for name, record in replay.get("input_files", {}).items():
        path = REPO_ROOT / name
        if path.is_file() and file_record(path) != normalize_record(record):
            errors.append(f"{study}: replay input differs from recorded provenance: {name}")
    if git_available():
        commit = replay.get("git_commit")
        if not commit:
            errors.append(f"{study}: replay source commit is missing")
        elif subprocess.run(
            ("git", "cat-file", "-e", f"{commit}^{{commit}}"), cwd=REPO_ROOT,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        ).returncode != 0:
            errors.append(f"{study}: replay source commit is unavailable in this checkout: {commit}")
    return errors


def analysis_warnings(study: str, manifest: Mapping[str, Any]) -> list[str]:
    restatement = manifest.get("restatement")
    if not isinstance(restatement, Mapping):
        return []
    recorded = restatement.get("analysis_pipeline_files", {})
    if not recorded:
        return [f"{study}: restatement provenance has no analysis pipeline hashes"]
    return [f"{study}: current analysis differs from last restatement: {msg}" for msg in _compare_record_maps(recorded, analysis_records(study))]


def begin_restatement(study: str, source_bundle: Path) -> RestatementContext:
    source_bundle = source_bundle.resolve()
    cfg = STUDIES[study]
    source_manifest = source_bundle.parent / cfg["manifest"]
    source_raw = source_bundle.parent / cfg["raw_csv"]
    if not source_manifest.is_file():
        raise ProvenanceError(f"restatement requires source manifest {source_manifest}")
    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    replay = manifest.get("replay_provenance")
    if not isinstance(replay, Mapping):
        raise ProvenanceError(
            "ERROR: source bundle predates immutable replay provenance. "
            "Use the explicit verified legacy migration or perform a full replay; "
            "statistical restatement cannot establish replay provenance."
        )
    errors = replay_errors(study, manifest)
    if errors:
        raise ProvenanceError(
            "ERROR: simulator/estimator implementation differs from replay provenance.\n"
            "A full replay is required; statistical restatement cannot update this bundle.\n"
            + "\n".join(errors)
        )
    return RestatementContext(
        study=study,
        source_bundle=source_bundle,
        source_manifest=source_manifest,
        source_raw_csv=source_raw,
        replay_provenance=copy.deepcopy(dict(replay)),
        source_bundle_sha256=sha256_file(source_bundle),
        source_manifest_sha256=sha256_file(source_manifest),
    )


def preserve_raw_rows(context: RestatementContext, destination: Path) -> None:
    if sha256_file(context.source_raw_csv) != context.replay_provenance["raw_rows_sha256"]:
        raise ProvenanceError("raw replay rows changed after restatement validation")
    if context.source_raw_csv.resolve() == destination.resolve():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(context.source_raw_csv, destination)


def finalize_restatement_manifest(
    study: str,
    manifest: dict[str, Any],
    context: RestatementContext,
) -> dict[str, Any]:
    result = dict(manifest)
    result["schema_version"] = SCHEMA_VERSION
    result["replay_provenance"] = copy.deepcopy(context.replay_provenance)
    result["restatement"] = {
        "git_commit": git_output("rev-parse", "HEAD"),
        "analysis_pipeline_files": analysis_records(study),
        "source_bundle_sha256": context.source_bundle_sha256,
        "source_manifest_sha256": context.source_manifest_sha256,
        "restated_at": utc_now(),
        "environment": environment_metadata(),
    }
    _replay_aliases(result)
    return result


def migrate_legacy_manifest(study: str, replay_commit: str) -> None:
    """One-time, auditable migration of a legacy *restated* committed bundle.

    This is not a replay and never claims to be one.  It succeeds only when Git
    history proves (a) the nominated historical manifest was a full generation,
    (b) every current replay dependency is unchanged since that commit, and (c)
    the committed normalized raw CSV is exactly the historical replay CSV with
    only the retired gate columns removed.
    """
    if not git_available():
        raise ProvenanceError("legacy replay migration requires full Git history")
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    manifest_path = out / cfg["manifest"]
    raw_path = out / cfg["raw_csv"]
    current = json.loads(manifest_path.read_text(encoding="utf-8"))
    if current.get("replay_provenance"):
        raise ProvenanceError(f"{study}: replay provenance is already migrated")

    manifest_repo_path = repo_name(manifest_path)
    raw_repo_path = repo_name(raw_path)
    historical = json.loads(git_bytes(replay_commit, manifest_repo_path).decode("utf-8"))
    command = [str(item) for item in historical.get("command", [])]
    if historical.get("git_commit") != replay_commit or "--restat-from" in command or historical.get("restated_from"):
        raise ProvenanceError(
            f"{study}: {replay_commit} does not contain a full-generation manifest"
        )

    paths = [repo_name(path) for path in implementation_paths(study)]
    changed = subprocess.run(
        ("git", "diff", "--quiet", replay_commit, "--", *paths), cwd=REPO_ROOT,
        check=False,
    ).returncode
    if changed != 0:
        raise ProvenanceError(
            "ERROR: simulator/estimator implementation differs from the nominated "
            "historical replay commit. A full replay is required."
        )

    historical_raw = git_bytes(replay_commit, raw_repo_path)
    if not normalized_raw_rows_equal_historical(raw_path, historical_raw):
        raise ProvenanceError(
            f"{study}: normalized raw rows are not the historical replay rows with "
            "only legacy gate columns removed"
        )

    old_analysis = {
        name: normalize_record(record)
        for name, record in current.get("analysis_pipeline_files", {}).items()
    }
    replay_environment = {
        "python": historical.get("python", "unrecorded"),
        "numpy": historical.get("numpy", "unrecorded"),
        "platform": historical.get("platform", "unrecorded"),
        "compiler": {
            "recorded": False,
            "note": "legacy replay predates compiler identity provenance",
        },
        "eigen": {
            "recorded": False,
            "note": "legacy replay predates explicit Eigen identity provenance",
        },
    }
    replay = {
        "git_commit": replay_commit,
        "implementation_files": implementation_records(study),
        "input_files": _input_records_from_legacy(historical),
        "raw_rows_sha256": sha256_file(raw_path),
        "raw_rows_size_bytes": raw_path.stat().st_size,
        "raw_rows_schema": "normalized statistical replay rows; deterministic simulator gate columns excluded",
        "workflow_run": None,
        "generated_at": None,
        "environment": replay_environment,
        "provenance_migration": {
            "kind": "verified_legacy_full_replay",
            "historical_manifest_commit": replay_commit,
            "historical_raw_rows_sha256": sha256_bytes(historical_raw),
            "normalization_only_columns_removed": list(GATE_FIELDS),
            "verified_at_git_commit": git_output("rev-parse", "HEAD"),
            "verified_at": utc_now(),
        },
    }
    current["schema_version"] = SCHEMA_VERSION
    current["replay_provenance"] = replay
    if old_analysis:
        current["restatement"] = {
            "git_commit": current.get("git_commit"),
            "analysis_pipeline_files": old_analysis,
            "source_bundle_sha256": sha256_file(out / cfg["json"]),
            "source_manifest_sha256": None,
            "restated_at": None,
            "environment": {
                "python": current.get("python", "unrecorded"),
                "numpy": current.get("numpy", "unrecorded"),
                "platform": current.get("platform", "unrecorded"),
            },
            "provenance_migration_note": "restatement metadata recovered from the legacy manifest; timestamp/source-manifest hash were not recorded",
        }
    _replay_aliases(current)
    manifest_path.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def initialize_fresh_manifest(study: str, manifest: dict[str, Any]) -> dict[str, Any]:
    result = dict(manifest)
    result["schema_version"] = SCHEMA_VERSION
    result["replay_provenance"] = create_replay_provenance_from_fresh_generation(study, result)
    result.pop("restatement", None)
    result["analysis_pipeline_files"] = _legacy_map(analysis_records(study))
    _replay_aliases(result)
    return result
