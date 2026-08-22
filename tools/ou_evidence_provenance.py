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

# Historical replay manifests hashed the entire multi-target OU-III Makefile.
# Removing unrelated diagnostic targets must not relabel the scientific replay,
# but changing anything that builds the simulator still must invalidate it.
# Keep the immutable historical record and special-case only this known legacy
# representation.  Every estimator/source dependency remains byte-for-byte.
_OU_III_MAKEFILE_REPO_NAME = "tests/kalman_ou_iii/Makefile"
_OU_III_HISTORICAL_MAKEFILE_SHA256 = (
    "d18838eb8c009408292fba8dcd09d94f97a13d01d743a9dce96ed30bb59b31eb"
)
_OU_III_SIM_BUILD_CONTRACT = """CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
CXXFLAGS = $(BASEFLAGS) -march=native
else
CXXFLAGS = $(BASEFLAGS)
endif
LDFLAGS  =
VPATH = $(REPO_ROOT)/src/util
kalman_ou_iii-sim: kalman_ou_iii-sim.o W3dSimCommon.o
$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.cpp
$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
"""
_BUILD_ASSIGNMENT_PREFIXES = (
    "CC =",
    "TEST_DIR :=",
    "REPO_ROOT :=",
    "EIGEN_DIR ?=",
    "EIGEN_CPPFLAGS :=",
    "BASEFLAGS =",
    "UNAME_S :=",
    "CXXFLAGS =",
    "LDFLAGS",
    "VPATH =",
)

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


def _ou_iii_sim_build_contract(text: str) -> str:
    """Project the multi-target Makefile onto replay-producing semantics."""
    lines = text.splitlines()
    selected: list[str] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        if (
            stripped.startswith(_BUILD_ASSIGNMENT_PREFIXES)
            or stripped.startswith("ifeq (")
            or stripped.startswith("else ifeq (")
            or stripped in {"else", "endif"}
        ):
            selected.append(stripped)
        if stripped.startswith("kalman_ou_iii-sim:"):
            selected.append(stripped)
            if index + 1 < len(lines) and lines[index + 1].startswith("\t"):
                selected.append(lines[index + 1].strip())
        if stripped == "%.o: %.cpp":
            selected.append(stripped)
            if index + 1 < len(lines) and lines[index + 1].startswith("\t"):
                selected.append(lines[index + 1].strip())
    return "\n".join(selected) + "\n"


def implementation_record_matches(
    name: str,
    expected: Mapping[str, Any],
    actual: Mapping[str, Any] | None = None,
) -> bool:
    """Compare replay implementation records without weakening source hashes.

    The only non-bytewise compatibility case is the historically over-broad
    OU-III Makefile record.  It is accepted only for the one known immutable
    historical hash and only while the simulator-producing projection remains
    exactly equal to the recorded build contract.  Fresh replay manifests still
    record the complete current Makefile through ``implementation_records``.
    """
    path = REPO_ROOT / name
    current = normalize_record(actual) if actual is not None else (
        file_record(path) if path.is_file() else None
    )
    if current is not None and normalize_record(expected) == current:
        return True
    if name != _OU_III_MAKEFILE_REPO_NAME:
        return False
    if normalize_record(expected)["sha256"] != _OU_III_HISTORICAL_MAKEFILE_SHA256:
        return False
    if not path.is_file():
        return False
    projected = _ou_iii_sim_build_contract(path.read_text(encoding="utf-8"))
    return projected == _OU_III_SIM_BUILD_CONTRACT


def implementation_record_errors(
    expected: Mapping[str, Mapping[str, Any]],
    actual: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[str]:
    """Compare replay dependencies using byte identity except the legacy Makefile."""
    actual_records = implementation_records("validation") if actual is None else actual
    errors: list[str] = []
    for name in sorted(set(expected) | set(actual_records)):
        if name not in expected:
            errors.append(f"new replay dependency not present in replay provenance: {name}")
        elif name not in actual_records:
            errors.append(f"recorded replay dependency missing from source tree: {name}")
        elif not implementation_record_matches(name, expected[name], actual_records[name]):
            errors.append(f"replay dependency differs from replay provenance: {name}")
    return errors


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


def _input_records_from_legacy(study: str, manifest: Mapping[str, Any]) -> dict[str, dict[str, object]]:
    # Legacy robustness manifests mixed replay inputs, simulator sources, and
    # Python analysis files under one source_files key. Preserve only actual
    # replay inputs here; implementation/build and analysis provenance have
    # their own explicit namespaces in schema v2.
    excluded = set(implementation_records(study)) | set(analysis_records(study))
    return {
        name: normalize_record(record)
        for name, record in manifest.get("source_files", {}).items()
        if name not in excluded
    }


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
        "input_files": _input_records_from_legacy(study, manifest),
        "raw_rows_sha256": sha256_file(raw),
        "raw_rows_size_bytes": raw.stat().st_size,
        "raw_rows_schema": "normalized statistical replay rows; deterministic simulator gate columns excluded",
        "workflow_run": workflow_metadata(),
        "generated_at": utc_now(),
        "environment": manifest.get("build_environment") or environment_metadata(),
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
        actual_impl = implementation_records(study)
        errors.extend(
            f"{study}: {msg}"
            for msg in implementation_record_errors(expected_impl, actual_impl)
        )
    for name, record in replay.get("input_files", {}).items():
        path = REPO_ROOT / name
        if path.is_file() and file_record(path) != normalize_record(record):
            errors.append(f"{study}: replay input differs from recorded provenance: {name}")
    commit = str(replay.get("git_commit", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        errors.append(f"{study}: replay source commit is missing or malformed")
    # A shallow Git checkout may not contain the historical replay object.
    # Hash validation above is authoritative for the committed source tree;
    # full Git history is required only for the explicit legacy migration.
    return errors


def analysis_warnings(study: str, manifest: Mapping[str, Any]) -> list[str]:
    restatement = manifest.get("restatement")
    if not isinstance(restatement, Mapping):
        return []
    recorded = restatement.get("analysis_pipeline_files", {})
    if not recorded:
        return [f"{study}: restatement provenance has no analysis pipeline hashes"]
    return [f"{study}: current analysis differs from last restatement: {msg}" for msg in _compare_record_maps(recorded, analysis_records(study))]


def _scalar_matches_csv(value: Any, csv_value: str) -> bool:
    if value is None:
        return csv_value == ""
    if isinstance(value, bool):
        return csv_value.lower() in ({"true", "1"} if value else {"false", "0"})
    if isinstance(value, (int, float)) and csv_value != "":
        try:
            lhs = float(value)
            rhs = float(csv_value)
            if lhs != lhs and rhs != rhs:  # NaN
                return True
            return lhs == rhs
        except ValueError:
            pass
    return str(value) == csv_value


def _assert_bundle_rows_match_raw_csv(source_bundle: Path, raw_csv: Path) -> None:
    bundle = json.loads(source_bundle.read_text(encoding="utf-8"))
    rows = bundle.get("raw_runs", [])
    with raw_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        csv_rows = list(reader)
    if len(rows) != len(csv_rows):
        raise ProvenanceError("bundle raw_runs count differs from canonical raw replay CSV")
    expected_fields = set(fields)
    for index, (row, csv_row) in enumerate(zip(rows, csv_rows)):
        extra = set(row) - expected_fields
        if extra:
            raise ProvenanceError(
                f"bundle raw_runs has fields outside canonical raw CSV at row {index}: {sorted(extra)}"
            )
        for field in fields:
            if field not in row:
                if csv_row[field] != "":
                    raise ProvenanceError(
                        f"bundle raw_runs omits non-empty canonical field at row {index}, field {field}"
                    )
                continue
            if not _scalar_matches_csv(row[field], csv_row[field]):
                raise ProvenanceError(
                    f"bundle raw_runs contradicts canonical raw replay CSV at row {index}, field {field}"
                )


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
    if not source_raw.is_file():
        errors.append(f"{study}: canonical raw replay CSV is missing: {source_raw}")
    if errors:
        raise ProvenanceError(
            "ERROR: simulator/estimator implementation differs from replay provenance.\n"
            "A full replay is required; statistical restatement cannot update this bundle.\n"
            + "\n".join(errors)
        )
    _assert_bundle_rows_match_raw_csv(source_bundle, source_raw)
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


def migrate_legacy_manifest(
    study: str,
    replay_commit: str,
    historical_evidence_commit: str,
) -> None:
    """One-time, auditable migration of a legacy *restated* committed bundle.

    `replay_commit` identifies the source tree that produced simulator outputs.
    `historical_evidence_commit` identifies the later commit that first stored
    the corresponding full-generation manifest and raw rows. Keeping these two
    identities separate is essential: an evidence-archive commit must never be
    misrepresented as the estimator revision that produced the replay.

    The migration succeeds only when Git history proves that:
      1. the historical evidence manifest says its replay source was
         `replay_commit` and is not itself a restatement;
      2. every replay-producing dependency in the current closure is unchanged
         since `replay_commit`; and
      3. the current normalized raw CSV is byte-identical to the historical
         full-replay CSV after removing only the retired gate columns.
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
    historical = json.loads(
        git_bytes(historical_evidence_commit, manifest_repo_path).decode("utf-8")
    )
    command = [str(item) for item in historical.get("command", [])]
    historical_restat = (
        historical.get("restated_from")
        or historical.get("protocol", {}).get("restated_from")
        or "--restat-from" in command
    )
    if historical.get("git_commit") != replay_commit or historical_restat:
        raise ProvenanceError(
            f"{study}: {historical_evidence_commit} is not a full-generation "
            f"bundle produced from replay source {replay_commit}"
        )

    paths = [repo_name(path) for path in implementation_paths(study)]
    changed = subprocess.run(
        ("git", "diff", "--quiet", replay_commit, "--", *paths),
        cwd=REPO_ROOT,
        check=False,
    ).returncode
    if changed != 0:
        changed_names = subprocess.run(
            ("git", "diff", "--name-only", replay_commit, "--", *paths),
            cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, check=False,
        ).stdout.strip()
        raise ProvenanceError(
            "ERROR: simulator/estimator implementation differs from replay "
            "provenance. A full replay is required; statistical migration "
            "cannot update this bundle."
            + (f"\nChanged replay dependencies:\n{changed_names}" if changed_names else "")
        )

    historical_raw = git_bytes(historical_evidence_commit, raw_repo_path)
    if not normalized_raw_rows_equal_historical(raw_path, historical_raw):
        raise ProvenanceError(
            f"{study}: normalized raw rows are not the historical full-replay "
            "rows with only legacy gate columns removed; a replay or forensic "
            "reconstruction is required"
        )

    old_analysis = {
        name: normalize_record(record)
        for name, record in current.get("analysis_pipeline_files", {}).items()
    }
    current_command = [str(item) for item in current.get("command", [])]
    current_is_restatement = bool(
        current.get("restated_from")
        or current.get("protocol", {}).get("restated_from")
        or "--restat-from" in current_command
    )
    historical_fields = csv.DictReader(
        io.StringIO(historical_raw.decode("utf-8"), newline="")
    ).fieldnames or []
    removed_gate_columns = [name for name in GATE_FIELDS if name in historical_fields]
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
        "input_files": _input_records_from_legacy(study, historical),
        "raw_rows_sha256": sha256_file(raw_path),
        "raw_rows_size_bytes": raw_path.stat().st_size,
        "raw_rows_schema": (
            "normalized statistical replay rows; deterministic simulator gate "
            "columns excluded"
        ),
        "workflow_run": None,
        "generated_at": None,
        "environment": replay_environment,
        "provenance_migration": {
            "kind": "verified_legacy_full_replay",
            "replay_source_commit": replay_commit,
            "historical_evidence_commit": historical_evidence_commit,
            "historical_raw_rows_sha256": sha256_bytes(historical_raw),
            "normalization_only_columns_removed": removed_gate_columns,
            "verified_at_git_commit": git_output("rev-parse", "HEAD"),
            "verified_at": utc_now(),
        },
    }
    current["schema_version"] = SCHEMA_VERSION
    current["replay_provenance"] = replay
    if old_analysis and current_is_restatement:
        current["restatement"] = {
            "git_commit": current.get("git_commit"),
            "analysis_pipeline_files": old_analysis,
            "source_bundle_sha256": None,
            "source_manifest_sha256": None,
            "restated_at": None,
            "environment": {
                "python": current.get("python", "unrecorded"),
                "numpy": current.get("numpy", "unrecorded"),
                "platform": current.get("platform", "unrecorded"),
            },
            "provenance_migration_note": (
                "restatement metadata recovered from the legacy manifest; "
                "source bundle/manifest hashes and timestamp were not recorded "
                "and are intentionally not fabricated"
            ),
        }
    else:
        current.pop("restatement", None)
    _replay_aliases(current)
    manifest_path.write_text(
        json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def initialize_fresh_manifest(study: str, manifest: dict[str, Any]) -> dict[str, Any]:
    result = dict(manifest)
    result["schema_version"] = SCHEMA_VERSION
    result["replay_provenance"] = create_replay_provenance_from_fresh_generation(study, result)
    result.pop("restatement", None)
    result["analysis_pipeline_files"] = _legacy_map(analysis_records(study))
    _replay_aliases(result)
    return result
