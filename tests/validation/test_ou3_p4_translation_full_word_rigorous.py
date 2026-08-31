from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_translation_full_word_rigorous as R


class P4RigorousReusedTranslationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = R.build()

    def test_validates(self):
        self.assertEqual(R.validate(self.d), [])

    def test_reuses_438_without_unsafe_setup(self):
        self.assertTrue(self.d["reused_from_PR_438"])
        self.assertFalse(self.d["unsafe_438_pre_interval_scale_setup_used"])
        for mode in ("H", "A"):
            row = self.d["modes"][mode]
            self.assertTrue(row["source_float_literals_rounded_as_binary32"])
            self.assertTrue(row["conditioned_scale_products_outward_rounded_from_first_operation"])
            self.assertFalse(row["pre_interval_binary64_scale_products_used"])

    def test_full_word_translation_is_strictly_wider(self):
        for mode in ("H", "A"):
            row = self.d["modes"][mode]
            self.assertGreater(row["complete_word_translation_margin_lower"], 0.0)
            self.assertGreater(row["margin_widening_factor_lower"], 1.0)

    def test_partial_result_stays_fail_closed(self):
        self.assertEqual(self.d["P4_USABLE_CERTIFICATE_STATUS"], "NOT_ESTABLISHED")


if __name__ == "__main__":
    unittest.main()
