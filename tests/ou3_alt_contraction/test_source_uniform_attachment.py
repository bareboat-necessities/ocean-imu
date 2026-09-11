"""ALT source-uniform attachment audit remains deliberately fail-closed."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'tools/stability')]
from tools.stability.ou3_alt_contraction import source_uniform_attachment as audit


class SourceUniformAttachmentTests(unittest.TestCase):
    def test_event_relation_executes_without_promoting_pair(self):
        d = audit.build()
        self.assertEqual(audit.validate(d), [])
        self.assertTrue(d['event_local_attachment_executes_without_startup_theorem'])
        self.assertTrue(d['inherited_theorem_builder']['startup_coupled'])
        self.assertFalse(d['paired_finite_increment_relation_closed'])
        self.assertFalse(d['ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED'])
        self.assertFalse(d['ALT_LIVE_PASS'])
        self.assertFalse(d['rank3_structure_safe_to_use']['full_state_reduction_permitted'])
        self.assertEqual(d['rank3_structure_safe_to_use']['thin_master_port_dimension_A21'], 33)
        self.assertEqual(d['rank3_structure_safe_to_use']['old_vec_increment_dimension_A21'], 105)
        missing = [k for k, v in d['paired_finite_increment_requirements'].items() if not v]
        self.assertIn('deltaN_and_deltaS_bound_to_same_paired_history', missing)
        self.assertIn('paired_guard_outcomes_or_hard_guard_graph', missing)
        self.assertIn('thin_product_ports_uN_uS_have_hard_same_history_product_graphs', missing)


if __name__ == '__main__':
    unittest.main()
