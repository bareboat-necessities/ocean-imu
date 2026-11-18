"""Bind startup replays to unchanged source bytes across compilation and scoring."""

from pathlib import Path
import subprocess

from sim_dataset import ROOT, digest


class SourceSnapshot:
    def __init__(self, root, paths):
        self.root = Path(root).resolve()
        self.hashes = {str(p.relative_to(self.root)): digest(p) for p in sorted(set(map(Path, paths)))}
        self.revision = self.current_revision()

    def current_revision(self):
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=self.root, text=True).strip()

    def verify(self):
        changed = [
            name
            for name, checksum in self.hashes.items()
            if not (self.root / name).is_file() or digest(self.root / name) != checksum
        ]
        if self.current_revision() != self.revision:
            changed.append("HEAD")
        if changed:
            raise ValueError("startup replay sources changed during execution: " + ", ".join(changed))


def startup_snapshot(script):
    # Bind all library dependencies, including transitive template headers.
    # Taking these hashes only after a replay could mislabel an older binary.
    sources = [p for p in (ROOT / "src").rglob("*") if p.suffix in (".h", ".hpp", ".cpp", ".tpp", ".c")]
    experiment = ROOT / "tests/startup_learning"
    sources.extend([*experiment.glob("*.cpp"), *experiment.glob("*.h"), experiment / "Makefile"])
    sources.extend(
        ROOT / name
        for name in (
            "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
            "tools/startup_learning.py",
            "tools/startup_handover_components.py",
            "tools/startup_live_retention.py",
            "tools/sim_dataset.py",
            "tools/startup_source_snapshot.py",
        )
    )
    sources.append(Path(script).resolve())
    return SourceSnapshot(ROOT, sources)
