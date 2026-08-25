import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_subdivision as S


class Ou3P5FirstAccelSubdivisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = S.build(source_pieces=2, cayley_axis_pieces=4)

    def test_subdivision_is_source_bound_and_validates(self):
        self.assertEqual(S.validate(self.d), [])
        self.assertTrue(self.d["source_generated_not_trajectory_fit"])
        self.assertFalse(self.d["source_replay_used"])
        self.assertFalse(self.d["filter_changed"])
        self.assertTrue(self.d["uses_v3_dependency_preserving_full_matrix_backend"])

    def test_splits_requested_dependencies_without_raising_correction_limit(self):
        self.assertEqual(self.d["source_parameter_subdivision"], ["tau", "sigma_aw", "R_S"])
        self.assertTrue(self.d["cayley_ball_direction_subdivided"])
        self.assertTrue(self.d["pseudo_due_and_not_due_branches_kept_separate"])
        self.assertEqual(self.d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(self.d["deployed_correction_limit_increased"])

    def test_diagnostic_never_promotes_whole_P5_word(self):
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertGreater(self.d["evaluated_child_count"], 0)
        if self.d["P5_FIRST_ACCEL_SOURCE_CAYLEY_SUBDIVISION"] != "PASS":
            self.assertIsNotNone(self.d["first_unclosed_child"])
            self.assertEqual(
                self.d["next_obligation"],
                "ACCEL_VECTOR_ORIENTATION_AND_STATE_DIRECTION_SUBDIVISION_REQUIRED",
            )


if __name__ == "__main__":
    unittest.main()
