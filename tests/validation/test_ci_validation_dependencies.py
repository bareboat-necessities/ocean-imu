"""Keep reusable validation jobs able to run the complete test suite.

The full-study fingerprint/reuse job has its own dependency installation. A
passing smoke or evidence-contract job therefore does not cover that path.
Only the standard library is needed to inspect the workflow's job blocks.
"""

from pathlib import Path
import re
import shlex
import unittest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/ou-validation.yml"
EVIDENCE_WORKFLOW = ROOT / ".github/workflows/evidence-contract.yml"
VALIDATION_COMMAND = re.compile(r"\bmake\s+-C\s+tests/validation\s+(?:evidence-)?test\b")
REQUIRED_PACKAGES = {
    "g++",
    "make",
    "libeigen3-dev",
    "python3-numpy",
    "python3-matplotlib",
    "python3-pandas",
    "python3-scipy",
}


def validation_job_prefixes(workflow):
    """Return each validation job's setup before its first test invocation."""
    jobs = workflow.split("\njobs:\n", 1)[1]
    blocks = re.split(r"(?m)^  ([A-Za-z_][A-Za-z0-9_-]*):\s*\n", jobs)
    prefixes = {}
    for job, body in zip(blocks[1::2], blocks[2::2]):
        command = VALIDATION_COMMAND.search(body)
        if command is not None:
            prefixes[job] = body[:command.start()]
    return prefixes


def installed_packages(prefix):
    """Read the actual apt install commands, not comments or later jobs."""
    commands = prefix.replace("\\\n", " ")
    packages = set()
    for line in commands.splitlines():
        line = line.split("#", 1)[0]
        install = re.search(r"\b(?:apt_get|apt-get)\b[^;\n]*?\binstall\b(?P<args>[^;\n]*)", line)
        if install is not None:
            packages.update(set(shlex.split(install.group("args"))) & REQUIRED_PACKAGES)
    return packages


class ValidationDependencyTests(unittest.TestCase):
    def test_every_validation_job_installs_the_full_test_dependencies(self):
        for workflow, expected in ((WORKFLOW, {"validate", "fingerprint", "commit"}),
                                   (EVIDENCE_WORKFLOW, {"evidence-contract"})):
            prefixes = validation_job_prefixes(workflow.read_text(encoding="utf-8"))
            self.assertEqual(set(prefixes), expected)
            for job, prefix in prefixes.items():
                with self.subTest(workflow=workflow.name, job=job):
                    self.assertFalse(
                        REQUIRED_PACKAGES - installed_packages(prefix),
                        f"{job} must install all validation build and import dependencies before running tests",
                    )

    def test_missing_compiler_or_eigen_is_not_covered_by_python_imports(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        for package in ("g++", "make", "libeigen3-dev"):
            with self.subTest(package=package):
                prefix = validation_job_prefixes(workflow)["commit"]
                # A removed prerequisite must be detected, even if all Python
                # packages (and the same package in another job) remain.
                prefix = prefix.replace(package, "")
                self.assertEqual(REQUIRED_PACKAGES - installed_packages(prefix), {package})

    def test_wrapped_install_command_includes_native_dependencies(self):
        prefix = "apt_get install -y g++ make libeigen3-dev " + "\\" + "\n  python3-numpy"
        self.assertEqual(installed_packages(prefix),
                         {"g++", "make", "libeigen3-dev", "python3-numpy"})

    def test_dependency_in_another_job_or_comment_does_not_cover_reuse(self):
        workflow = """
jobs:
  validate:
    run: apt_get install python3-pandas libeigen3-dev g++
  fingerprint:
    # apt_get install python3-pandas libeigen3-dev g++ is not an executed command.
    run: apt_get install python3-numpy
    test: make -C tests/validation evidence-test
    late: apt_get install python3-pandas libeigen3-dev g++
"""
        prefixes = validation_job_prefixes(workflow)
        self.assertEqual(set(prefixes), {"fingerprint"})
        self.assertEqual(installed_packages(prefixes["fingerprint"]), {"python3-numpy"})


if __name__ == "__main__":
    unittest.main()
