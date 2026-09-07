from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_physical_sea3_membership as M  # noqa: E402


class PhysicalSea3MembershipTests(unittest.TestCase):
    def test_third_harmonic_fails_even_largest_declared_response(self):
        d = M.particle_response_test("12/5", "4", "6/5", "2")
        self.assertFalse(d["response_envelope_satisfied"])
        self.assertEqual(d["particle_response_norm_squared_exact"], "2")
        self.assertEqual(d["maximum_allowed_norm_squared_exact"], "1")
        self.assertEqual(d["squared_gap_exact"], "1")

    def test_fundamental_is_not_falsely_excluded_by_that_test(self):
        d = M.particle_response_test("4/5", "4", "6/5", "2")
        self.assertTrue(d["response_envelope_satisfied"])
        self.assertEqual(d["maximum_allowed_norm_squared_exact"], "16")

    def test_steeper_rolloff_cannot_repair_high_frequency_violation(self):
        d = M.particle_response_test("12/5", "4", "6/5", "3")
        self.assertFalse(d["response_envelope_satisfied"])
        self.assertEqual(d["maximum_allowed_norm_squared_exact"], "1/4")

    def test_invalid_response_parameters_are_rejected(self):
        for args in (("0", "4", "6/5", "2"), ("12/5", "-1", "6/5", "2"),
                     ("12/5", "4", "6/5", "1")):
            with self.assertRaises(ValueError):
                M.particle_response_test(*args)


if __name__ == "__main__":
    unittest.main()
