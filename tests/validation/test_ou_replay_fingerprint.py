import importlib.util
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "ou_replay_fingerprint.py"
SPEC = importlib.util.spec_from_file_location("ou_replay_fingerprint", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
fingerprint = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fingerprint)


class ReplayFingerprintClassificationTests(unittest.TestCase):
    def test_expected_source_and_workflow_types_are_included(self):
        names = [
            "src/a.c",
            "src/a.cc",
            "src/a.cpp",
            "src/a.cxx",
            "src/a.h",
            "src/a.hh",
            "src/a.hpp",
            "src/a.hxx",
            "src/a.ino",
            "src/a.S",
            "src/a.asm",
            "src/a.ipp",
            "src/a.tpp",
            "src/a.inc",
            "tools/a.py",
            "tools/a.sh",
            "tools/a.bash",
            "tools/a.zsh",
            "tools/a.pl",
            "tools/a.rb",
            "tools/a.lua",
            "tools/a.js",
            "tools/a.ts",
            ".github/workflows/build.yml",
            "config/build.yaml",
            "cmake/toolchain.cmake",
            "rules/common.mk",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(fingerprint.is_replay_source("100644", name))

    def test_expected_build_files_are_included(self):
        names = [
            "Makefile",
            "tests/foo/Makefile",
            "CMakeLists.txt",
            "Dockerfile",
            "meson.build",
            "meson_options.txt",
            "BUILD",
            "BUILD.bazel",
            "WORKSPACE",
            "WORKSPACE.bazel",
            "MODULE.bazel",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(fingerprint.is_replay_source("100644", name))

    def test_git_executable_mode_is_always_included(self):
        self.assertTrue(
            fingerprint.is_replay_source("100755", "tools/no-recognized-extension")
        )

    def test_generated_evidence_and_prose_do_not_self_invalidate(self):
        names = [
            "reports/results/ou_validation/ou_validation.json",
            "reports/results/ou_robustness/ou_robustness_raw.csv",
            "doc/kalman_ou_iii/w3d-ou-validation-results-generated.tex-part",
            "README.md",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertFalse(fingerprint.is_replay_source("100644", name))


if __name__ == "__main__":
    unittest.main()
