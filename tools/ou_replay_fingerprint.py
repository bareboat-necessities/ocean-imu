#!/usr/bin/env python3
"""Content fingerprint for deciding whether OU simulator evidence must be replayed.

The gate is deliberately conservative. It hashes every tracked file under
``tests/`` regardless of extension, every tracked source/script or build/workflow
file elsewhere, every tracked file whose Git mode is executable, and the exact
downloaded simulation-data ZIP.

This is intentionally broader than the simulator dependency closure. False
positive full replays are acceptable; reusing evidence after a potentially
replay-affecting repository change is not.

The fingerprint is independent of the enclosing Git commit SHA and of generated
evidence. A documentation/evidence-only commit outside ``tests/`` therefore does
not invalidate a scientifically identical replay. Files under ``tests/`` are
always included because test configuration and study parameters may use
arbitrary file extensions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 1
ALGORITHM = "sha256-framed-v1"

# Deliberately broad. Do not replace this with include/reachability analysis:
# the replay gate is designed to fail toward unnecessary recomputation rather
# than toward accidental reuse of stale evidence.
SOURCE_SUFFIXES = {
    # C/C++/Objective-C/Arduino/assembly and source fragments.
    ".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx",
    ".m", ".mm", ".ino", ".s", ".asm", ".ipp", ".tpp", ".inc",
    # Scripts and executable-language sources.
    ".py", ".sh", ".bash", ".zsh", ".fish", ".pl", ".rb", ".lua",
    ".js", ".mjs", ".cjs", ".ts", ".mts", ".cts", ".java", ".kt",
    ".kts", ".rs", ".go", ".swift",
    # Workflow/build configuration that can change how executable code runs.
    ".yml", ".yaml", ".mk", ".cmake", ".ac", ".am",
}
SOURCE_BASENAMES = {
    "Makefile",
    "CMakeLists.txt",
    "Dockerfile",
    "meson.build",
    "meson_options.txt",
    "BUILD",
    "BUILD.bazel",
    "WORKSPACE",
    "WORKSPACE.bazel",
    "MODULE.bazel",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_ls_files() -> list[tuple[str, str]]:
    completed = subprocess.run(
        ("git", "ls-files", "-s", "-z"),
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "git ls-files failed: "
            + completed.stderr.decode("utf-8", errors="replace").strip()
        )

    entries: list[tuple[str, str]] = []
    for raw in completed.stdout.split(b"\0"):
        if not raw:
            continue
        meta, path_bytes = raw.split(b"\t", 1)
        mode = meta.split(b" ", 1)[0].decode("ascii")
        path = path_bytes.decode("utf-8", errors="surrogateescape")
        entries.append((mode, path))
    return entries


def is_replay_source(mode: str, relative_path: str) -> bool:
    path = Path(relative_path)
    return (
        relative_path.startswith("tests/")
        or path.suffix.lower() in SOURCE_SUFFIXES
        or path.name in SOURCE_BASENAMES
        or mode == "100755"
    )


def replay_source_entries() -> list[tuple[str, str]]:
    return sorted(
        (entry for entry in _git_ls_files() if is_replay_source(*entry)),
        key=lambda item: item[1],
    )


def _frame(digest, *parts: bytes) -> None:
    for part in parts:
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)


def compute_fingerprint(simulation_zip: Path) -> dict[str, object]:
    simulation_zip = simulation_zip.resolve()
    if not simulation_zip.is_file():
        raise FileNotFoundError(simulation_zip)

    digest = hashlib.sha256()
    _frame(digest, b"ocean-imu-ou-replay-fingerprint", str(SCHEMA_VERSION).encode())

    files: list[dict[str, object]] = []
    for mode, relative_path in replay_source_entries():
        path = REPO_ROOT / relative_path
        if not path.is_file():
            raise FileNotFoundError(path)
        file_sha = sha256_file(path)
        files.append({"path": relative_path, "mode": mode, "sha256": file_sha})
        _frame(
            digest,
            b"repo-file",
            mode.encode("ascii"),
            relative_path.encode("utf-8", errors="surrogateescape"),
            bytes.fromhex(file_sha),
        )

    zip_sha = sha256_file(simulation_zip)
    _frame(digest, b"simulation-data-zip", bytes.fromhex(zip_sha))

    return {
        "schema_version": SCHEMA_VERSION,
        "algorithm": ALGORITHM,
        "fingerprint": digest.hexdigest(),
        "simulation_data": {
            "file": simulation_zip.name,
            "sha256": zip_sha,
        },
        "tracked_replay_file_count": len(files),
        "tracked_replay_files": files,
    }


def _comparison_view(record: dict[str, object]) -> tuple[object, ...]:
    simulation = record.get("simulation_data")
    simulation_sha = simulation.get("sha256") if isinstance(simulation, dict) else None
    return (
        record.get("schema_version"),
        record.get("algorithm"),
        record.get("fingerprint"),
        simulation_sha,
    )


def write_record(path: Path, record: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def check_record(path: Path, current: dict[str, object]) -> bool:
    if not path.is_file():
        print(f"replay fingerprint missing: {path}", file=sys.stderr)
        return False
    try:
        recorded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"cannot read replay fingerprint {path}: {exc}", file=sys.stderr)
        return False
    if not isinstance(recorded, dict):
        print(f"replay fingerprint is not a JSON object: {path}", file=sys.stderr)
        return False
    if _comparison_view(recorded) != _comparison_view(current):
        print(
            "OU replay fingerprint changed: "
            f"recorded={recorded.get('fingerprint')} current={current.get('fingerprint')}",
            file=sys.stderr,
        )
        return False
    print(f"OU replay fingerprint unchanged: {current['fingerprint']}")
    return True


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--simulation-zip",
        type=Path,
        required=True,
        help="exact sim-data-files.zip used by the studies",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write", type=Path, help="write the current fingerprint JSON")
    action.add_argument("--check", type=Path, help="return 0 only when this record matches")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    current = compute_fingerprint(args.simulation_zip)
    if args.write is not None:
        write_record(args.write, current)
        print(f"wrote OU replay fingerprint {current['fingerprint']} to {args.write}")
        return 0
    if args.check is not None:
        return 0 if check_record(args.check, current) else 1
    json.dump(current, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
