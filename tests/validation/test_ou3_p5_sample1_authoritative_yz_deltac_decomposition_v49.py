import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_authoritative_yz_deltac_decomposition_v49 as V49


class Sample1AuthoritativeYzDeltaCDecompositionV49Tests(unittest.TestCase):
    def test_authoritative_witness_and_schema_are_frozen(self):
        self.assertEqual(V49.WITNESS, (0, 0, 23))
        self.assertEqual(V49.SCHEMA, 4900)

    def test_deltac_terms_are_outward_and_sum_to_candidate(self):
        d = V49._deltac_terms(dP=2.0, dH=3.0, htheta=4.0, row_norm=5.0)
        self.assertGreaterEqual(d["projected_DeltaP_Htheta_upper"], 8.0)
        self.assertGreaterEqual(d["nominal_Ptheta_row_DeltaH_upper"], 15.0)
        self.assertGreaterEqual(d["mixed_DeltaP_DeltaH_upper"], 6.0)
        self.assertGreaterEqual(d["theta_aw_cross_block_parent_upper"], 2.0)
        self.assertGreaterEqual(d["row_DeltaC_candidate_upper"], 31.0)
        self.assertIn(d["dominant_term"], d)
        self.assertTrue(math.isfinite(d["dominant_fraction"]))
        self.assertGreaterEqual(d["dominant_fraction"], 0.0)
        self.assertLessEqual(d["dominant_fraction"], 1.0)

    def test_invalid_terms_fail_closed(self):
        with self.assertRaises(ValueError):
            V49._deltac_terms(dP=-1.0, dH=1.0, htheta=1.0, row_norm=1.0)


if __name__ == "__main__":
    unittest.main()
