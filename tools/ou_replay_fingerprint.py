#!/usr/bin/env python3
"""Content fingerprints for deciding whether OU evidence must be replayed.

The replay gate is deliberately conservative. It hashes every tracked file under
``tests/`` regardless of extension, every tracked source/script or build/workflow
file elsewhere, every tracked file whose Git mode is executable, and the exact
downloaded simulation-data ZIP.

A second fingerprint hashes every scientific result below ``reports/results/``
except ``reports/results/readme/``.  That directory is a derived presentation
mirror populated from successful build artifacts; its bytes are not primary
validation evidence and must not make an unchanged study look tampered with.
The record that stores the two fingerprints lives outside the results tree, so
the scientific digest has no self-reference exception.

This is intentionally broader than the simulator dependency closure. False
positive full replays are acceptable; reusing evidence after a potentially
replay-affecting repository or scientific-evidence change is not.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = REPO_ROOT / "reports" / "results"
RESULTS_FINGERPRINT_EXCLUDED_TOP_LEVEL = frozenset({"readme"})
SCHEMA_VERSION = 3
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
    ".yml", ".yaml", ".mk", ".make", ".cmake", ".ninja", ".ac", ".am",
}
SOURCE_BASENAMES = {
    "Makefile",
    "CMakeLists.txt",
    "CMakePresets.json",
    "CMakeUserPresets.json",
    "Dockerfile",
    "meson.build",
    "meson_options.txt",
    "BUILD",
    "BUILD.bazel",
    "WORKSPACE",
    "WORKSPACE.bazel",
    "MODULE.bazel",
    "pyproject.toml",
    "setup.cfg",
    "tox.ini",
    "platformio.ini",
    "Pipfile",
    "Pipfile.lock",
    "poetry.lock",
    "package.json",
    "package-lock.json",
    "yarn.lock",
    "pnpm-lock.yaml",
    "Cargo.toml",
    "Cargo.lock",
    "go.mod",
    "go.sum",
    "build.gradle",
    "settings.gradle",
    "gradle.properties",
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


def _is_build_filename(name: str) -> bool:
    lower = name.lower()
    return (
        name in SOURCE_BASENAMES
        or name.startswith("Makefile")
        or name.startswith("Dockerfile")
        or (lower.startswith("requirements") and lower.endswith(".txt"))
    )


def is_replay_source(mode: str, relative_path: str) -> bool:
    path = Path(relative_path)
    return (
        relative_path.startswith("tests/")
        or relative_path.startswith(".github/")
        or path.suffix.lower() in SOURCE_SUFFIXES
        or _is_build_filename(path.name)
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


def is_results_fingerprint_excluded(path: Path, results_root: Path = RESULTS_ROOT) -> bool:
    """Return whether ``path`` belongs to a derived presentation-only subtree."""

    relative = path.relative_to(results_root)
    return bool(relative.parts) and relative.parts[0] in RESULTS_FINGERPRINT_EXCLUDED_TOP_LEVEL


def compute_results_fingerprint(results_root: Path = RESULTS_ROOT) -> dict[str, object]:
    """Hash protected result files/symlinks and their paths.

    ``readme/`` is intentionally excluded: it mirrors successful build output
    for repository presentation and is not scientific evidence.
    """

    results_root = results_root.resolve()
    if not results_root.is_dir():
        raise FileNotFoundError(results_root)

    digest = hashlib.sha256()
    _frame(digest, b"ocean-imu-ou-results-fingerprint", str(SCHEMA_VERSION).encode())

    files: list[dict[str, object]] = []
    for path in sorted(results_root.rglob("*"), key=lambda item: item.as_posix()):
        if is_results_fingerprint_excluded(path, results_root):
            continue
        if path.is_symlink():
            target = os.readlink(path)
            content_sha = hashlib.sha256(target.encode("utf-8")).hexdigest()
            kind = "symlink"
        elif path.is_file():
            content_sha = sha256_file(path)
            kind = "file"
        else:
            continue

        try:
            relative_path = path.relative_to(REPO_ROOT).as_posix()
        except ValueError:
            relative_path = path.relative_to(results_root).as_posix()

        files.append(
            {"path": relative_path, "kind": kind, "sha256": content_sha}
        )
        _frame(
            digest,
            b"results-entry",
            kind.encode("ascii"),
            relative_path.encode("utf-8", errors="surrogateescape"),
            bytes.fromhex(content_sha),
        )

    return {
        "root": "reports/results" if results_root == RESULTS_ROOT.resolve() else str(results_root),
        "fingerprint": digest.hexdigest(),
        "file_count": len(files),
        "files": files,
    }


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
    results = compute_results_fingerprint()

    return {
        "schema_version": SCHEMA_VERSION,
        "algorithm": ALGORITHM,
        "replay_fingerprint": digest.hexdigest(),
        "simulation_data": {
            "file": simulation_zip.name,
            "sha256": zip_sha,
        },
        "tracked_replay_file_count": len(files),
        "tracked_replay_files": files,
        "results_fingerprint": results["fingerprint"],
        "results_file_count": results["file_count"],
        "results_files": results["files"],
    }


def _comparison_view(record: dict[str, object]) -> tuple[object, ...]:
    simulation = record.get("simulation_data")
    simulation_sha = simulation.get("sha256") if isinstance(simulation, dict) else None
    return (
        record.get("schema_version"),
        record.get("algorithm"),
        record.get("replay_fingerprint"),
        simulation_sha,
        record.get("results_fingerprint"),
    )


def write_record(path: Path, record: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def check_record(path: Path, current: dict[str, object]) -> bool:
    if not path.is_file():
        print(f"OU evidence fingerprint missing: {path}", file=sys.stderr)
        return False
    try:
        recorded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"cannot read OU evidence fingerprint {path}: {exc}", file=sys.stderr)
        return False
    if not isinstance(recorded, dict):
        print(f"OU evidence fingerprint is not a JSON object: {path}", file=sys.stderr)
        return False
    if _comparison_view(recorded) != _comparison_view(current):
        print(
            "OU replay/results fingerprint changed: "
            f"replay recorded={recorded.get('replay_fingerprint')} "
            f"current={current.get('replay_fingerprint')}; "
            f"results recorded={recorded.get('results_fingerprint')} "
            f"current={current.get('results_fingerprint')}",
            file=sys.stderr,
        )
        return False
    print(
        "OU replay/results fingerprints unchanged: "
        f"replay={current['replay_fingerprint']} "
        f"results={current['results_fingerprint']}"
    )
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
        print(
            "wrote OU evidence fingerprints "
            f"replay={current['replay_fingerprint']} "
            f"results={current['results_fingerprint']} to {args.write}"
        )
        return 0
    if args.check is not None:
        return 0 if check_record(args.check, current) else 1
    json.dump(current, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
