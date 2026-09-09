from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1] / "kalman_ou_iii"
TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_p4_uniform_closure as C


class UniformClosureTests(unittest.TestCase):
    def test_contract_closes_without_covariance_entry_assumption(self):
        d = C.build()
        self.assertEqual(C.validate(d), [])
        self.assertTrue(d["P4_MOTION_PASS"])
        self.assertTrue(d["P4_PASS"])
        self.assertFalse(d["shipping_covariance_used_as_entry_membership_test"])
        self.assertFalse(d["independent_P_H_R_K_boxes_used"])
        self.assertTrue(d["BIAS1_SOURCE_ADMISSION_PASS"])
        self.assertTrue(d["finite_precision_enclosure_closed"])
        for mode in ("H18", "A21"):
            m = d["modes"][mode]
            self.assertGreater(m["certified_hard_entry_scale"], 0.0)
            self.assertTrue(m["consecutive_compatible_storage_closed"])
            self.assertFalse(m["finite_precision_zero_floor_contraction_claimed"])


if __name__ == "__main__":
    unittest.main()
