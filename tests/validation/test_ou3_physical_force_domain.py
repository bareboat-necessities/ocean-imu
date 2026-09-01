from pathlib import Path
import copy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_physical_force_domain as F


class PhysicalForceDomainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = F.build()

    def test_marine_force_domain_is_derived_and_valid(self):
        self.assertEqual(F.validate(self.d), [])
        self.assertAlmostEqual(self.d["gravity_mps2"], 9.80665, places=12)
        self.assertAlmostEqual(self.d["non_gravitational_cog_acceleration_norm_upper_mps2"], 4.0, places=12)
        self.assertAlmostEqual(self.d["non_gravitational_cog_acceleration_fraction_g"], 0.4078864851911713, places=15)
        self.assertAlmostEqual(self.d["derived_specific_force_norm_lower_mps2"], 5.80665, places=12)
        self.assertAlmostEqual(self.d["derived_specific_force_norm_upper_mps2"], 13.80665, places=12)
        self.assertFalse(self.d["old_independent_30_mps2_cap_used"])
        self.assertFalse(self.d["impact_slam_in_scope"])
        self.assertFalse(self.d["lever_arm_enabled"])

    def test_total_specific_force_is_not_wave_acceleration(self):
        self.assertEqual(self.d["primary_physical_assumption"], "||a_non-grav,CoG|| <= a_max")
        self.assertLess(self.d["derived_specific_force_norm_lower_fraction_g"], 1.0)
        self.assertGreater(self.d["derived_specific_force_norm_upper_fraction_g"], 1.0)
        self.assertLess(self.d["derived_specific_force_norm_upper_fraction_g"], 1.5)

    def test_validator_rejects_old_three_g_envelope(self):
        d = copy.deepcopy(self.d)
        d["old_independent_30_mps2_cap_used"] = True
        self.assertIn("old_independent_30_mps2_cap_used is not false", F.validate(d))


if __name__ == "__main__":
    unittest.main()
