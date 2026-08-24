import json
import sys
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_explicit_information_word_certificate as CERT


class Ou3ExplicitInformationWordCertificateTests(unittest.TestCase):
    def test_source_uniform_H_A_information_margins_are_strict(self):
        d = CERT.build()
        self.assertEqual(CERT.validate(d), [])
        self.assertFalse(d["source_generated_not_trajectory_fit"] is False)
        self.assertEqual(d["continuous_linear_information_certificate"], "PASS")
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "LINEAR_ONLY")
        for mode in ("H", "A"):
            row = d["modes"][mode]
            self.assertGreater(row["Sigma_lambda_min_lower"], 0.0)
            self.assertGreater(row["Sigma_lambda_max_upper"], row["Sigma_lambda_min_lower"])
            self.assertGreater(row["word_noise_Omega_lambda_min_lower"], 0.0)
            self.assertGreater(row["relative_Riccati_injection_margin_lower"], 0.0)
            self.assertLess(row["relative_Riccati_injection_margin_lower"], 1.0)
            self.assertEqual(row["prefix_information_gain_upper"], 1.0)

    def test_certificate_depends_on_declared_physical_upper_bounds(self):
        base = json.loads(CERT.DEFAULT_DOMAIN.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "domain.json"
            broken = json.loads(json.dumps(base))
            broken["normal_live"].pop("specific_force_norm_upper_mps2")
            p.write_text(json.dumps(broken), encoding="utf-8")
            with self.assertRaises(KeyError):
                CERT.build(p)


if __name__ == "__main__":
    unittest.main()
