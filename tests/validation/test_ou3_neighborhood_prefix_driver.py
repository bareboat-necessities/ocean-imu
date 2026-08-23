import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_neighborhood_prefix_driver", ROOT / "tools" / "ou3_neighborhood_prefix_driver.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class PrefixDriverTests(unittest.TestCase):
    def test_prefix_keeps_header_and_first_sample_after_stop(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            src = root / "source.csv"
            dst = root / "out.csv"
            src.write_text(
                "time_s,x\n0.000,1\n0.005,2\n0.010,3\n0.015,4\n0.020,5\n"
            )
            meta = mod.write_prefix(src, dst, 0.011)
            lines = dst.read_text().splitlines()
            self.assertEqual(lines, ["time_s,x", "0.000,1", "0.005,2", "0.010,3", "0.015,4"])
            self.assertAlmostEqual(meta["last_copied_time_s"], 0.015)

    def test_prefix_rejects_too_short_source(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            src = root / "source.csv"
            dst = root / "out.csv"
            src.write_text("time_s,x\n0.000,1\n0.005,2\n")
            with self.assertRaises(RuntimeError):
                mod.write_prefix(src, dst, 1.0)


if __name__ == "__main__":
    unittest.main()
