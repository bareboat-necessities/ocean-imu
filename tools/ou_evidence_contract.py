#!/usr/bin/env python3
"""Validate and normalize committed OU statistical evidence.

The contract distinguishes two provenance layers:

* replay_provenance: immutable source commit, replay-producing implementation
  dependency closure, replay inputs, and the canonical normalized raw-row hash;
* restatement: later analysis/editorial context used to regenerate derived
  statistics, tables, and publication text from the already-existing rows.

Only a genuine full simulator regeneration may create replay provenance. A
statistical restatement may never replace it or make rows from an older
estimator appear current.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping

import ou_evidence_provenance as provenance

REPO_ROOT = provenance.REPO_ROOT
STUDIES = provenance.STUDIES
GATE_FIELDS = set(provenance.GATE_FIELDS)

# Compatibility exports used by existing tests and tooling.
sha256_file = provenance.sha256_file
file_record = provenance.file_record
repo_name = provenance.repo_name
implementation_closure = provenance.implementation_closure

# Historical replay manifests hashed the entire multi-target OU-III test
# Makefile.  That made unrelated test-target edits appear to invalidate the
# simulator replay.  Preserve the immutable replay record, but accept the known
# historical Makefile record when the current replay-producing build semantics
# are exactly the same.  Changes to compiler flags, include paths, simulator
# dependencies, link command, or object compilation command still invalidate
# the replay.
_OU_III_MAKEFILE = REPO_ROOT / "tests" / "kalman_ou_iii" / "Makefile"
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


def _historical_makefile_mismatch_is_auxiliary_only(
    study: str,
    manifest: Mapping[str, Any],
    error: str,
) -> bool:
    expected_error = (
        f"{study}: replay dependency differs from replay provenance: "
        f"{_OU_III_MAKEFILE_REPO_NAME}"
    )
    if error != expected_error:
        return False
    replay = manifest.get("replay_provenance")
    if not isinstance(replay, Mapping):
        return False
    record = replay.get("implementation_files", {}).get(_OU_III_MAKEFILE_REPO_NAME)
    if not isinstance(record, Mapping):
        return False
    if provenance.normalize_record(record)["sha256"] != _OU_III_HISTORICAL_MAKEFILE_SHA256:
        return False
    if not _OU_III_MAKEFILE.is_file():
        return False
    current = _ou_iii_sim_build_contract(_OU_III_MAKEFILE.read_text(encoding="utf-8"))
    return current == _OU_III_SIM_BUILD_CONTRACT


def _replay_errors(study: str, manifest: Mapping[str, Any]) -> list[str]:
    errors = provenance.replay_errors(study, manifest)
    return [
        error
        for error in errors
        if not _historical_makefile_mismatch_is_auxiliary_only(study, manifest, error)
    ]


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _manifest(study: str) -> tuple[Path, dict[str, Any]]:
    cfg = STUDIES[study]
    path = Path(cfg["dir"]) / cfg["manifest"]
    return path, json.loads(path.read_text(encoding="utf-8"))


def _refresh_result_files(study: str, manifest: dict[str, Any]) -> None:
    out = Path(STUDIES[study]["dir"])
    refreshed: dict[str, dict[str, object]] = {}
    for name in manifest.get("result_files", {}):
        path = out / name
        if path.is_file():
            refreshed[name] = {
                "sha256": provenance.sha256_file(path),
                "bytes": path.stat().st_size,
            }
    manifest["result_files"] = refreshed


def _statistical_schema_errors(study: str) -> list[str]:
    cfg = STUDIES[study]
    out = Path(cfg["dir"])
    errors: list[str] = []

    raw = out / cfg["raw_csv"]
    with raw.open(newline="", encoding="utf-8") as handle:
        fields = csv.DictReader(handle).fieldnames or []
    leaked = GATE_FIELDS.intersection(fields)
    if leaked:
        errors.append(f"{study}: legacy gate columns present: {sorted(leaked)}")

    bundle = json.loads((out / cfg["json"]).read_text(encoding="utf-8"))
    if any(GATE_FIELDS.intersection(row) for row in bundle.get("raw_runs", [])):
        errors.append(f"{study}: legacy gate fields present in raw_runs")
    protocol = bundle.get("protocol", {})
    if protocol.get("simulator_regression_gates_exported") is not False:
        errors.append(f"{study}: statistical evidence gate policy is not explicit")
    if "all completed replays" not in protocol.get("replay_inclusion_rule", ""):
        errors.append(f"{study}: replay inclusion rule is missing or ambiguous")

    macro_name = cfg["macro"]
    if macro_name:
        text = (out / macro_name).read_text(encoding="utf-8")
        if "OUValidationGatePasses" in text or "OUValidationGateFailures" in text:
            errors.append("validation: legacy gate-count macros present")
    return errors


def _result_inventory_errors(study: str, manifest: Mapping[str, Any]) -> list[str]:
    out = Path(STUDIES[study]["dir"])
    errors: list[str] = []
    for name, record in manifest.get("result_files", {}).items():
        path = out / name
        if not path.is_file():
            errors.append(f"{study}: committed result missing: {name}")
            continue
        expected = provenance.normalize_record(record)
        actual = provenance.file_record(path)
        if expected != actual:
            errors.append(f"{study}: committed result hash mismatch: {name}")
    return errors


def check(study: str, *, require_current_analysis: bool = False) -> list[str]:
    """Return hard contract violations.

    Replay validity is always enforced. A later source edit to analysis code is
    not allowed to invalidate the scientific replay itself; callers that need
    byte-for-byte reproducibility of the latest restatement can additionally
    request current-analysis matching.
    """
    errors = _statistical_schema_errors(study)
    _, manifest = _manifest(study)
    if manifest.get("schema_version") != provenance.SCHEMA_VERSION:
        errors.append(f"{study}: evidence manifest schema_version must be {provenance.SCHEMA_VERSION}")
    errors.extend(_replay_errors(study, manifest))
    errors.extend(_result_inventory_errors(study, manifest))
    if require_current_analysis:
        errors.extend(provenance.analysis_warnings(study, manifest))
    return errors


def normalize(study: str) -> None:
    """Normalize statistical schema without ever mutating replay provenance."""
    provenance.normalize_statistical_schema(study)
    path, manifest = _manifest(study)
    if manifest.get("replay_provenance"):
        _refresh_result_files(study, manifest)
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return

    # This path is only valid immediately after a full regeneration at HEAD.
    manifest = provenance.initialize_fresh_manifest(study, manifest)
    _refresh_result_files(study, manifest)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _auto(studies: tuple[str, ...]) -> int:
    # In an unpacked source archive, --auto becomes a read-only check. There is
    # no Git identity with which to create new replay provenance, and there
    # should not be: archives validate what was committed rather than inventing
    # a source commit.
    if not provenance.git_available():
        errors = [error for study in studies for error in check(study)]
        if errors:
            print("\n".join(errors))
            return 1
        print("OU evidence contract is current (archive/hash mode; .git unavailable)")
        return 0

    for study in studies:
        _, manifest = _manifest(study)
        if manifest.get("replay_provenance"):
            continue
        command = [str(item) for item in manifest.get("command", [])]
        restated = manifest.get("restated_from") or manifest.get("protocol", {}).get("restated_from")
        if restated or "--restat-from" in command:
            print(
                f"{study}: legacy restated bundle has no immutable replay provenance; "
                "--auto will not convert it. Use explicit verified migration or a full replay."
            )
            return 1

    # Only manifests that provenance.initialize_fresh_manifest can prove are
    # freshly generated at the current HEAD are normalized/stamped.
    for study in studies:
        _, manifest = _manifest(study)
        if not manifest.get("replay_provenance"):
            try:
                normalize(study)
            except provenance.ProvenanceError as exc:
                print(str(exc))
                return 1

    errors = [error for study in studies for error in check(study)]
    if errors:
        print("\n".join(errors))
        return 1
    print("OU evidence contract is current")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("validation", "robustness", "all"), default="all")
    parser.add_argument(
        "--require-current-analysis",
        action="store_true",
        help="also require current analysis scripts to match the most recent restatement provenance",
    )
    parser.add_argument(
        "--replay-commit",
        help="source commit whose simulator/estimator implementation produced the replay rows; migration only",
    )
    parser.add_argument(
        "--historical-evidence-commit",
        help="commit containing the original full-generation manifest/raw rows for --migrate-legacy",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="read-only verification")
    mode.add_argument(
        "--auto",
        action="store_true",
        help="initialize provenance only for a freshly generated full replay at HEAD; otherwise check",
    )
    mode.add_argument(
        "--migrate-legacy",
        action="store_true",
        help="one-time verified migration of legacy restated manifests; requires full Git history",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    studies = tuple(STUDIES) if args.study == "all" else (args.study,)

    if args.migrate_legacy:
        if not args.replay_commit or not args.historical_evidence_commit:
            raise SystemExit(
                "--migrate-legacy requires both --replay-commit and "
                "--historical-evidence-commit"
            )
        try:
            for study in studies:
                provenance.migrate_legacy_manifest(
                    study,
                    replay_commit=args.replay_commit,
                    historical_evidence_commit=args.historical_evidence_commit,
                )
        except provenance.ProvenanceError as exc:
            print(str(exc))
            return 1
        errors = [error for study in studies for error in check(study)]
        if errors:
            print("\n".join(errors))
            return 1
        print("Migrated legacy OU evidence provenance without replaying simulator rows")
        return 0

    if args.auto:
        return _auto(studies)

    errors = [
        error
        for study in studies
        for error in check(study, require_current_analysis=args.require_current_analysis)
    ]
    if errors:
        print("\n".join(errors))
        return 1

    for study in studies:
        for warning in provenance.analysis_warnings(study, _manifest(study)[1]):
            if not args.require_current_analysis:
                print(f"NOTE: {warning}")
    mode = "archive/hash mode" if not provenance.git_available() else "Git-aware mode"
    print(f"OU evidence contract is current ({mode})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
