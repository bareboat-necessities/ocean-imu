"""Keep the CMake route parseable and in step with its sources.

The Make-based suites never read CMakeLists.txt, so a corrupted file (a
literal "\\n" once joined two add_executable commands) could sit unnoticed
while every Make target passed. The full configure+build runs in the build
workflow's cmake job; this check configures locally whenever cmake exists.
"""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
CMAKE = ROOT / "CMakeLists.txt"
BUILD_WORKFLOW = ROOT / ".github" / "workflows" / "build.yml"


class CMakeBuildContract(unittest.TestCase):
    def test_no_escaped_newline_corruption(self):
        for number, line in enumerate(CMAKE.read_text(encoding="utf-8").splitlines(), 1):
            with self.subTest(line=number):
                self.assertNotIn("\\n", line)
                self.assertLessEqual(line.count("add_executable("), 1, line)

    def test_every_listed_source_exists(self):
        text = CMAKE.read_text(encoding="utf-8")
        paths = re.findall(r'"\$\{(OCEAN_IMU_TESTS_DIR|OCEAN_IMU_SRC_DIR)\}/([^"]+)"', text)
        self.assertTrue(paths)
        for variable, relative in paths:
            base = ROOT / ("tests" if variable == "OCEAN_IMU_TESTS_DIR" else "src")
            with self.subTest(path=relative):
                self.assertTrue((base / relative).is_file())

    def test_build_workflow_configures_and_builds_with_cmake(self):
        workflow = BUILD_WORKFLOW.read_text(encoding="utf-8")
        start = workflow.index("\n  cmake:\n")
        job = workflow[start:]
        self.assertIn("cmake -S . -B build", job)
        self.assertIn("cmake --build build", job)
        self.assertNotIn("continue-on-error", job)

    @unittest.skipUnless(shutil.which("cmake"), "cmake is not installed")
    def test_configure_succeeds(self):
        with tempfile.TemporaryDirectory() as build:
            result = subprocess.run(
                ["cmake", "-S", str(ROOT), "-B", build],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
