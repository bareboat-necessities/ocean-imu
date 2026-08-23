import importlib.util
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou3_source_noise_certificate.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("ou3_source_noise_certificate", TOOL)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class Ou3SourceNoiseCertificateTests(unittest.TestCase):
    def test_source_generated_standardized_gaussian_bound(self):
        tool = load_tool()
        report = tool.build_certificate()
        self.assertTrue(report["source_generated_not_trajectory_fit"])
        z = report["standardized_increment"]
        self.assertEqual(z["dimension"], 18)
        self.assertTrue(z["covariance_upper_identity"])
        self.assertEqual(z["Sigma_bar_norm_upper"], 1.0)
        self.assertEqual(z["s2_upper"], 18.0)
        self.assertEqual(z["s4_upper"], 360.0)

    def test_physical_scales_match_deployed_validation_noise(self):
        tool = load_tool()
        p = tool.build_certificate()["physical_scales"]
        self.assertAlmostEqual(p["imu_dt_s"], 0.005, places=12)
        self.assertAlmostEqual(p["mag_odr_hz"], 25.0, places=12)
        self.assertAlmostEqual(p["acc_white_std_mps2"], 1.51e-3 * 9.80665, places=12)
        self.assertAlmostEqual(p["gyro_white_std_radps"], 0.00157, places=12)
        self.assertAlmostEqual(p["acc_bias_rw_increment_std_mps2"], 0.0005 * math.sqrt(0.005), places=12)
        self.assertAlmostEqual(p["gyro_bias_rw_increment_std_radps"], 1e-5 * math.sqrt(0.005), places=12)
        self.assertAlmostEqual(p["mag_white_std_uT"], 0.80, places=12)
        self.assertAlmostEqual(p["mag_bias_rw_increment_std_uT"], 0.01 * math.sqrt(1.0 / 25.0), places=12)


if __name__ == "__main__":
    unittest.main()
