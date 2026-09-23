"""Keep the native device-rate regression tied to the deployed sketch policy."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]


class StationaryDeviceContract(unittest.TestCase):
    def test_deployed_pseudo_cadence(self):
        for family in ("ou2", "ou3", "tfg"):
            directory = f"atomS3R_ins_kalman_{family}" if family != "tfg" else "atomS3R_ins_tfg"
            text = (ROOT / "sensors/full_marine_ins" / directory / f"{directory}.ino").read_text()
            # OU sketches use the fixed 15 ms cadence; TFG keeps the tau-scaled one.
            call = ("ff.setTauScaledPseudoUpdateCadence(false);" if family != "tfg"
                    else "fusion_.setTauScaledPseudoCadence(true);")
            with self.subTest(family=family):
                self.assertEqual(text.count(call), 1)
                # The call must be live code, not swallowed by a comment.
                line = next(l for l in text.splitlines() if call in l)
                self.assertTrue(line.strip().startswith(call), line)

    def test_shared_regression_rules_preserve_native_build_flags(self):
        text = (ROOT / "tests/common/StationaryDeviceRegression.mk").read_text()
        self.assertIn("stationary_device-test: stationary_device-test.o", text)
        self.assertIn("$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)", text)

    def test_all_three_native_suites_run_the_regression(self):
        for family in ("kalman_ou_ii", "kalman_ou_iii", "kalman_tfg"):
            path = ROOT / "tests" / family
            with self.subTest(family=family):
                self.assertIn("../common/StationaryDeviceRegression.mk stationary_device-test",
                              (path / "run_tests.sh").read_text())
                self.assertIn("./stationary_device-test", (path / "run_tests.sh").read_text())


if __name__ == "__main__":
    unittest.main()
