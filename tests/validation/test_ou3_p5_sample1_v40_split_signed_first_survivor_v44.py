import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_v40_split_signed_first_survivor_v44 as V44


class Sample1V40SplitSignedFirstSurvivorV44Tests(unittest.TestCase):
    def test_intersection_is_componentwise_and_nonworsening(self):
        a = [Interval(-2.0, 2.0), Interval(-1.0, 3.0), Interval(4.0, 8.0)]
        b = [Interval(-1.0, 1.0), Interval(0.0, 4.0), Interval(5.0, 7.0)]
        c = V44._intersect_boxes(a, b)
        self.assertIsNotNone(c)
        self.assertEqual([(x.lo, x.hi) for x in c], [(-1.0, 1.0), (0.0, 3.0), (5.0, 7.0)])

    def test_disjoint_component_makes_source_cell_incompatible(self):
        a = [Interval(-1.0, 0.0), Interval(-1.0, 1.0), Interval(-1.0, 1.0)]
        b = [Interval(0.5, 1.0), Interval(-1.0, 1.0), Interval(-1.0, 1.0)]
        self.assertIsNone(V44._intersect_boxes(a, b))

    def test_witness_and_targets_are_frozen(self):
        self.assertEqual(V44.WITNESS, (0, 0, 23))
        self.assertAlmostEqual(V44.V41_REFERENCE_FIRST_Q, 8.344528951460543, places=12)
        self.assertEqual(V44.Q_TARGET, 8.0)

    def test_contract_is_fail_closed(self):
        self.assertEqual(V44.SCHEMA, 4400)
        self.assertEqual(V44.Q_TARGET, 8.0)


if __name__ == "__main__":
    unittest.main()
