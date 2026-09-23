"""Keep the native device-rate regression tied to the deployed sketch policy."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]


class StationaryDeviceContract(unittest.TestCase):
    def test_deployed_pseudo_cadence(self):
        for family in ("ou2", "ou3", "tfg"):
            directory = f"atomS3R_ins_kalman_{family}" if family != "tfg" else "atomS3R_ins_tfg"
            text = (ROOT / "sensors/full_marine_ins" / directory / f"{directory}.ino").read_text()
            call = ("ff.setTauScaledPseudoUpdateCadence(false);" if family != "tfg"
                    else "fusion_.setTauScaledPseudoCadence(false);")
            with self.subTest(family=family):
                self.assertEqual(text.count(call), 1)

    def test_all_three_native_suites_run_the_regression(self):
        for family in ("kalman_ou_ii", "kalman_ou_iii", "kalman_tfg"):
            path = ROOT / "tests" / family
            with self.subTest(family=family):
                self.assertIn("stationary_device-test", (path / "Makefile").read_text())
                self.assertIn("./stationary_device-test", (path / "run_tests.sh").read_text())


if __name__ == "__main__":
    unittest.main()
