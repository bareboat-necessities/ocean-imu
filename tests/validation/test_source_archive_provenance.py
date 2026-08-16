"""Evidence-contract behavior in GitHub/source archives without a .git directory."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


class SourceArchiveProvenanceTests(unittest.TestCase):
    def _copy_archive_fixture(self, destination: Path) -> None:
        # Copy exactly the material the evidence contract needs.  Deliberately
        # do not copy .git: this models GitHub source ZIPs, release archives,
        # Zenodo artifacts, and ordinary copied source trees.
        for relative in (
            "src",
            "tests/kalman_ou_ii",
            "tests/kalman_ou_iii",
            "reports/results/ou_validation",
            "reports/results/ou_robustness",
        ):
            shutil.copytree(REPO_ROOT / relative, destination / relative)
        (destination / "tools").mkdir(parents=True)
        for name in (
            "ou_evidence_contract.py",
            "ou_evidence_provenance.py",
            "ou_validation.py",
            "ou_robustness.py",
        ):
            shutil.copy2(REPO_ROOT / "tools" / name, destination / "tools" / name)

    def _run(self, root: Path, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            (sys.executable, *args),
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_committed_evidence_verifies_without_git_metadata(self):
        with tempfile.TemporaryDirectory(prefix="ocean-imu-source-archive-") as temporary:
            root = Path(temporary)
            self._copy_archive_fixture(root)
            result = self._run(root, "tools/ou_evidence_contract.py", "--check")
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("archive/hash mode", result.stdout)

    def test_changed_estimator_cannot_be_laundered_by_restat_then_auto_in_archive(self):
        with tempfile.TemporaryDirectory(prefix="ocean-imu-source-archive-adversarial-") as temporary:
            root = Path(temporary)
            self._copy_archive_fixture(root)
            header = root / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
            header.write_text(
                header.read_text(encoding="utf-8") + "\n// adversarial provenance mutation\n",
                encoding="utf-8",
            )

            restat = self._run(
                root,
                "tools/ou_validation.py",
                "--restat-from",
                "reports/results/ou_validation/ou_validation.json",
                "--output-dir",
                "reports/results/ou_validation",
            )
            self.assertNotEqual(restat.returncode, 0, restat.stdout)
            self.assertIn("full replay is required", restat.stdout.lower())

            auto = self._run(root, "tools/ou_evidence_contract.py", "--auto")
            self.assertNotEqual(auto.returncode, 0, auto.stdout)
            self.assertIn("replay dependency differs", auto.stdout.lower())


if __name__ == "__main__":
    unittest.main()
