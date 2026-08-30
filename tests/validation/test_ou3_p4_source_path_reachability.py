from pathlib import Path
import sys, unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_source_path_reachability as R

class P4SourcePathReachabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.d=R.build()

    def test_validates(self):
        self.assertEqual(R.validate(self.d),[])

    def test_source_only_and_not_promoted(self):
        self.assertTrue(self.d['source_only'])
        self.assertFalse(self.d['trajectory_replay_used'])
        self.assertFalse(self.d['usable_P4_promoted'])

    def test_graph_is_nonempty(self):
        self.assertGreater(self.d['partition']['states'],0)
        self.assertGreater(self.d['transition_edges'],0)
        self.assertGreater(self.d['strongly_connected_components'],0)

    def test_old_worst_corner_is_explicit(self):
        self.assertGreater(self.d['old_worst_corner_state_count'],0)
        self.assertGreaterEqual(self.d['old_worst_corner_recurrent_state_count'],0)

if __name__=='__main__': unittest.main()
