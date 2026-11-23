"""Compare the firmware's serial writes to the original snprintf wire formats."""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class NmeaOutputTest(unittest.TestCase):
    def test_sentence_bytes(self):
        here = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / "nmea"
            (Path(tmp) / "Arduino.h").touch()  # The test provides the Serial sink.
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17", "-O2",
                            "-Wall", "-Wextra", "-Werror", "-I", tmp, "-I", str(here.parents[1] / "src"),
                            str(here / "nmea_output_test.cpp"), "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)


if __name__ == "__main__":
    unittest.main()
