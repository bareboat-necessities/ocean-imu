#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import symmetric_positive_definite_ldlt
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_segment as S
import ou3_p3_frozen_full_matrix_translation as F


class CorrelatedTranslationSegmentTests(unittest.TestCase):
    def test_contract_is_bound_to_frozen_p2_interface(self):
        d = S.build(representative_nodes=(137,), representative_gaps=(13,))
        self.assertEqual(S.validate(d), [])
        self.assertTrue(d["P2_correlation_interface_consumed"])
        self.assertEqual(d["P2_correlation_interface_version"], CORR.INTERFACE_VERSION)
        self.assertTrue(d["tau_sigma_R_S_same_source_cell_per_segment"])
        self.assertFalse(d["sigma_dependent_state_rescaling_used"])
        self.assertFalse(d["P3_PROMOTED"])

    def test_zero_start_segment_is_strict_full_matrix(self):
        rt = CORR.runtime()
        rows = S.segment_images(F._mat_zero(), 137, 13, rt, x_subcells=2)
        self.assertGreater(len(rows), 0)
        for row in rows:
            P = row["posterior"]
            self.assertEqual(len(P), 4)
            self.assertTrue(all(len(r) == 4 for r in P))
            self.assertTrue(symmetric_positive_definite_ldlt(P)[0])

    def test_legal_pair_shift_uses_staged_source_for_following_segment(self):
        rt = CORR.runtime()
        c = next(i for i, out in enumerate(rt["union_successors"]) if out)
        s = min(rt["union_successors"][c])
        gap = next(g for g in rt["gaps"] if CORR.successors(s, g, rt))
        t = CORR.successors(s, gap, rt)[0]
        tr = CORR.transition(c, s, gap, t, rt)
        self.assertEqual(tr["following_segment_applied_node"], s)
        rows = S.segment_images(F._mat_zero(), s, gap, rt, x_subcells=1)
        self.assertGreater(len(rows), 0)
        self.assertTrue(all(symmetric_positive_definite_ldlt(r["posterior"])[0] for r in rows))


if __name__ == "__main__":
    unittest.main()
