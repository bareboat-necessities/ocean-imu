"""Keep reusable validation jobs able to import the complete test suite.

The full-study fingerprint/reuse job has its own dependency installation. A
passing smoke or evidence-contract job therefore does not cover that path.
Only the standard library is needed to inspect the workflow's job blocks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/ou-validation.yml"
VALIDATION_COMMAND = re.compile(r"\bmake\s+-C\s+tests/validation\s+(?:evidence-)?test\b")
REQUIRED_PACKAGES = {
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


def installed_python_packages(prefix):
    """Read the actual apt install commands, not comments or later jobs."""
    commands = prefix.replace("\\\n", " ")
    packages = set()
    for line in commands.splitlines():
        line = line.split("#", 1)[0]
        if re.search(r"\b(?:apt_get|apt-get)\s+install\b", line):
            packages.update(re.findall(r"\bpython3-[a-z0-9-]+\b", line))
    return packages


class ValidationDependencyTests(unittest.TestCase):
    def test_every_validation_job_installs_the_full_python_test_dependencies(self):
        prefixes = validation_job_prefixes(WORKFLOW.read_text(encoding="utf-8"))
        self.assertEqual(set(prefixes), {"validate", "fingerprint", "commit"})
        for job, prefix in prefixes.items():
            with self.subTest(job=job):
                self.assertFalse(
                    REQUIRED_PACKAGES - installed_python_packages(prefix),
                    f"{job} must install all validation imports before running tests",
                )

    def test_dependency_in_another_job_or_comment_does_not_cover_reuse(self):
        workflow = """
jobs:
  validate:
    run: apt_get install python3-pandas
  fingerprint:
    # apt_get install python3-pandas is not an executed command.
    run: apt_get install python3-numpy
    test: make -C tests/validation evidence-test
    late: apt_get install python3-pandas
"""
        prefixes = validation_job_prefixes(workflow)
        self.assertEqual(set(prefixes), {"fingerprint"})
        self.assertEqual(installed_python_packages(prefixes["fingerprint"]), {"python3-numpy"})


if __name__ == "__main__":
    unittest.main()
