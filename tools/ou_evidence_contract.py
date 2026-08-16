#!/usr/bin/env python3
"""Validate and normalize committed OU statistical evidence.

The contract distinguishes two provenance layers:

* replay_provenance: immutable source commit, replay-producing implementation
  dependency closure, replay inputs, and the canonical normalized raw-row hash;
* restatement: later analysis/editorial context used to regenerate derived
  statistics, tables, and publication text from the already-existing rows.

Only a genuine full simulator regeneration may create replay provenance.  A
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

    Replay validity is always enforced.  A later source edit to analysis code is
    not allowed to invalidate the scientific replay itself; callers that need
    byte-for-byte reproducibility of the latest restatement can additionally
    request current-analysis matching.
    """
    errors = _statistical_schema_errors(study)
    _, manifest = _manifest(study)
    if manifest.get("schema_version") != provenance.SCHEMA_VERSION:
        errors.append(f"{study}: evidence manifest schema_version must be {provenance.SCHEMA_VERSION}")
    errors.extend(provenance.replay_errors(study, manifest))
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
        # Replay aliases are preserved by the provenance module; do not stamp
        # current implementation hashes over them here.
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return

    # This path is only valid immediately after a full regeneration at HEAD.
    manifest = provenance.initialize_fresh_manifest(study, manifest)
    _refresh_result_files(study, manifest)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _auto(studies: tuple[str, ...]) -> int:
    # In an unpacked source archive, --auto becomes a read-only check.  There is
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
        help="also require the current analysis scripts to match the most recent restatement provenance",
    )
    parser.add_argument(
        "--replay-commit",
        help="historical full-replay commit used only with --migrate-legacy",
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
        if not args.replay_commit:
            raise SystemExit("--migrate-legacy requires --replay-commit")
        try:
            for study in studies:
                provenance.migrate_legacy_manifest(study, args.replay_commit)
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
