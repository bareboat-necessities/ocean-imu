import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_v41_authoritative_split_signed_v45 as V45


class Sample1V41AuthoritativeSplitSignedV45Tests(unittest.TestCase):
    def test_archived_v41_witness_is_frozen(self):
        self.assertEqual(V45.WITNESS, (0, 0, 23))
        self.assertAlmostEqual(V45.V41_Q_CURRENT, 0.6415212986499801, places=14)
        self.assertAlmostEqual(V45.V41_Q_POST, 8.344528951460543, places=12)
        self.assertEqual(V45.Q_TARGET, 8.0)

    def test_reference_match_is_strict(self):
        self.assertTrue(V45._matches(V45.V41_Q_CURRENT, V45.V41_Q_CURRENT))
        self.assertFalse(V45._matches(V45.V41_Q_CURRENT + 1.0e-8, V45.V41_Q_CURRENT))

    def test_schema_is_distinct(self):
        self.assertEqual(V45.SCHEMA, 4500)


if __name__ == "__main__":
    unittest.main()
