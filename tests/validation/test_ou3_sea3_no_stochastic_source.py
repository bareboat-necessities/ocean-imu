from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_complete_source as source  # noqa: E402


class Sea3NoStochasticSourceTest(unittest.TestCase):
    def test_stochastic_calculation_is_forcing_only(self):
        d = source.build()
        self.assertEqual(source.validate(d), [])
        s = d["stochastic_forcing_corollary"]
        self.assertFalse(s["used_to_generate_P3_source_words"])
        self.assertFalse(s["used_to_prune_homogeneous_P3_family"])
        self.assertTrue(s["configured_Racc_Rmag_remain_in_every_covariance_update"])
        self.assertFalse(d["no_fallback_generators"]["gaussian_good_event_source_generator"])
        self.assertFalse(d["no_fallback_generators"]["spectral_moment_only_source_generator"])
        self.assertFalse(d["no_fallback_generators"]["arbitrary_bounded_input_source_generator"])


if __name__ == "__main__":
    unittest.main()
