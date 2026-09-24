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

    def test_tfg_device_tests_are_separate_executed_commands(self):
        directory = ROOT / "tests/kalman_tfg"
        commands = [line.strip() for line in (directory / "run_tests.sh").read_text().splitlines()]
        for name in ("tfg_device_heave-test", "tfg_device_tempcomp-test"):
            with self.subTest(name=name):
                self.assertEqual(commands.count("./" + name), 1)
                self.assertIn(f"{name}: {name}.o", (directory / "Makefile").read_text())

    def test_stationary_north_uses_each_familys_deployed_cadence(self):
        source = (ROOT / "tests/common/StationaryDeviceRegression.h").read_text()
        tfg = source.split("using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;", 1)[1]
        tfg, ou = tfg.split("#else", 1)
        self.assertIn("deployed_fixed_cadence = false;", tfg)
        self.assertIn("deployed_fixed_cadence = true;", ou)
        self.assertIn("replay(deployed_fixed_cadence,0.03f,0.025f,0.0f,false)", source)

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
