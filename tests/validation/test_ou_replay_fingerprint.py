import importlib.util
from pathlib import Path
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "ou_replay_fingerprint.py"
SPEC = importlib.util.spec_from_file_location("ou_replay_fingerprint", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
fingerprint = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fingerprint)


class ReplayFingerprintClassificationTests(unittest.TestCase):
    def test_every_file_under_tests_is_included_regardless_of_extension(self):
        names = [
            "tests/kalman_ou_iii/parameters.txt",
            "tests/kalman_ou_iii/fixture.csv",
            "tests/validation/config.json",
            "tests/validation/reference.dat",
            "tests/validation/README.md",
            "tests/custom/no-extension",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(fingerprint.is_replay_source("100644", name))

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
            ".github/dependabot/config.json",
            "config/build.yaml",
            "cmake/toolchain.cmake",
            "rules/common.mk",
            "rules/config.make",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(fingerprint.is_replay_source("100644", name))

    def test_expected_build_files_and_variants_are_included(self):
        names = [
            "Makefile",
            "Makefile.local",
            "Makefile.release",
            "CMakeLists.txt",
            "CMakePresets.json",
            "Dockerfile",
            "Dockerfile.ci",
            "meson.build",
            "meson_options.txt",
            "BUILD",
            "BUILD.bazel",
            "WORKSPACE",
            "WORKSPACE.bazel",
            "MODULE.bazel",
            "pyproject.toml",
            "setup.cfg",
            "platformio.ini",
            "requirements.txt",
            "requirements-ci.txt",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(fingerprint.is_replay_source("100644", name))

    def test_git_executable_mode_is_always_included(self):
        self.assertTrue(
            fingerprint.is_replay_source("100755", "tools/no-recognized-extension")
        )

    def test_generated_evidence_and_prose_outside_tests_do_not_self_invalidate_replay(self):
        names = [
            "reports/results/ou_validation/ou_validation.json",
            "reports/results/ou_robustness/ou_robustness_raw.csv",
            "reports/ou_evidence_fingerprint.json",
            "doc/kalman_ou_iii/w3d-ou-validation-results-generated.tex-part",
            "README.md",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertFalse(fingerprint.is_replay_source("100644", name))

    def test_results_fingerprint_covers_every_file_and_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "results"
            nested = root / "nested"
            nested.mkdir(parents=True)
            (root / "a.json").write_text('{"value": 1}\n', encoding="utf-8")
            (nested / "fixture.dat").write_bytes(b"abc")

            first = fingerprint.compute_results_fingerprint(root)
            self.assertEqual(first["file_count"], 2)

            (nested / "fixture.dat").write_bytes(b"abcd")
            second = fingerprint.compute_results_fingerprint(root)
            self.assertNotEqual(first["fingerprint"], second["fingerprint"])

            (nested / "new-extensionless-file").write_bytes(b"x")
            third = fingerprint.compute_results_fingerprint(root)
            self.assertNotEqual(second["fingerprint"], third["fingerprint"])
            self.assertEqual(third["file_count"], 3)


if __name__ == "__main__":
    unittest.main()
