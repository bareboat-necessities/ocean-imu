from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_physical_admissibility as physical  # noqa: E402
import ou3_brmm_complete_source as source  # noqa: E402


class BrmmCompactDomainTest(unittest.TestCase):
    def test_compactness_is_theorem_domain_property(self):
        p = physical.build()
        self.assertEqual(physical.validate(p), [])
        self.assertTrue(p["reference_parameter_domain_compact"])
        self.assertTrue(p["reference_compact_transition_relation"])
        self.assertTrue(p["P3_may_not_replace_compact_BRMM_with_independent_bounds"])

        d = source.build()
        self.assertEqual(source.validate(d), [])
        sea = d["reference_surface_family"]
        self.assertTrue(sea["parameter_domain_compact"])
        self.assertTrue(sea["reference_compact_transition_relation"])
        self.assertTrue(sea["independent_H_T_extrema_forbidden"])
        self.assertTrue(sea["independent_partition_height_maxima_forbidden"])


if __name__ == "__main__":
    unittest.main()
