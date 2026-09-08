"""Keep the canonical physical hypothesis consistent across repository sources."""
from pathlib import Path
import re
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[2]


class BrmmNamingTest(unittest.TestCase):
    def test_retired_source_identifier_is_absent_from_tracked_text_and_paths(self):
        retired = re.compile("s" + "ea" + "3", re.IGNORECASE)
        paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
        failures = []
        for name in filter(None, paths):
            p = ROOT / name
            if not p.is_file():
                continue
            if retired.search(name):
                failures.append(name)
            try:
                data = p.read_bytes()
                if b"\0" in data:
                    continue
                content = data.decode("utf-8")
            except UnicodeDecodeError:
                continue
            if retired.search(content):
                failures.append(name)
        self.assertEqual([], failures)


if __name__ == "__main__":
    unittest.main()
