#!/usr/bin/env python3
"""Reject ephemeral provenance and noncanonical JSON in OU-III result artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


FORBIDDEN_KEYS = {
    "generated_at",
    "generated_at_utc",
    "timestamp",
    "timestamp_utc",
    "build_id",
    "run_id",
    "run_number",
    "run_attempt",
    "workflow_id",
    "workflow_run_id",
    "artifact_id",
    "ci_run_id",
}


def walk_keys(value, path: str = "$"):
    if isinstance(value, dict):
        for key, child in value.items():
            key_text = str(key)
            yield path, key_text
            yield from walk_keys(child, f"{path}.{key_text}")
    elif isinstance(value, list):
        for i, child in enumerate(value):
            yield from walk_keys(child, f"{path}[{i}]")


def validate_json(path: Path) -> list[str]:
    failures: list[str] = []
    text = path.read_text(encoding="utf-8")
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        return [f"{path}: invalid JSON: {exc}"]

    for parent, key in walk_keys(payload):
        if key.lower() in FORBIDDEN_KEYS:
            failures.append(f"{path}: ephemeral key {parent}.{key}")

    canonical = json.dumps(payload, indent=2, sort_keys=True)
    if text != canonical:
        failures.append(
            f"{path}: JSON bytes are not canonical json.dumps(indent=2, sort_keys=True) output"
        )
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("result_dir", type=Path)
    args = ap.parse_args()
    root = args.result_dir.resolve()
    paths = sorted(root.rglob("*.json"))
    if not paths:
        raise SystemExit(f"no JSON result files found under {root}")

    failures: list[str] = []
    for path in paths:
        failures.extend(validate_json(path))
    if failures:
        print("OU-III result JSON determinism contract: FAIL")
        for failure in failures:
            print(f"- {failure}")
        return 2
    print(f"OU-III result JSON determinism contract: PASS ({len(paths)} JSON files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
