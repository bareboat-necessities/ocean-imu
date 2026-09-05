from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_physical_admissibility as phys  # noqa: E402


class Sea3PhysicalAdmissibilityTest(unittest.TestCase):
    def test_dnv_peak_steepness_piecewise_law(self) -> None:
        self.assertEqual(phys.peak_steepness_limit(1.0), 1.0 / 15.0)
        self.assertEqual(phys.peak_steepness_limit(8.0), 1.0 / 15.0)
        self.assertEqual(phys.peak_steepness_limit(15.0), 1.0 / 25.0)
        self.assertEqual(phys.peak_steepness_limit(30.0), 1.0 / 25.0)
        expected_mid = 0.5 * (1.0 / 15.0 + 1.0 / 25.0)
        self.assertTrue(math.isclose(phys.peak_steepness_limit(11.5), expected_mid))
        self.assertTrue(
            math.isclose(
                phys.peak_steepness_limit(math.nextafter(8.0, math.inf)),
                1.0 / 15.0,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                phys.peak_steepness_limit(math.nextafter(15.0, -math.inf)),
                1.0 / 25.0,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        )

    def test_height_period_coupling_blocks_unphysical_rectangular_corners(self) -> None:
        d = phys.build()
        self.assertEqual(phys.validate(d), [])
        g = d["gravity_mps2"]

        h1 = phys.significant_height_limit_from_peak_steepness(1.0, g)
        h8 = phys.significant_height_limit_from_peak_steepness(8.0, g)
        self.assertLess(h1, 0.11)
        self.assertGreater(h8, 6.6)
        self.assertLess(h8, 6.7)
        self.assertFalse(phys.partition_admissible(8.5, 1.0, g))
        self.assertTrue(phys.partition_admissible(0.1, 1.0, g))

        c = d["three_partition_contract"]
        self.assertTrue(c["independent_H_r_and_T_p_rectangular_extrema_forbidden"])
        self.assertTrue(c["independent_three_partition_H_maxima_forbidden"])
        self.assertEqual(c["total_Hs_upper_m"], 8.5)

    def test_parameter_subcertificate_retains_compact_sea3_without_faking_realization(self) -> None:
        d = phys.build()
        self.assertEqual(phys.validate(d), [])
        self.assertFalse(d["SEA0_full_certificate_promoted"])
        self.assertTrue(d["SEA3_parameter_domain_compact"])
        self.assertTrue(d["compact_transition_relation_is_theorem_domain"])
        self.assertTrue(d["this_subcertificate_refines_but_does_not_rectangularize_SEA3"])
        self.assertTrue(d["P3_may_not_replace_compact_SEA3_with_independent_bounds"])
        self.assertFalse(d["finite_window_realization_enclosed"])
        self.assertFalse(d["left_language_inclusion_closed"])
        self.assertTrue(
            d["external_basis"]["multimodal_partition_extension_is_our_conservative_theorem_choice"]
        )


if __name__ == "__main__":
    unittest.main()
