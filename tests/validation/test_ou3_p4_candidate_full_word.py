#!/usr/bin/env python3
import json
import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p4_candidate_full_word as C
import ou3_p4_p5_entrance_search_domain as E
import ou3_p4_operation_matched_sector_certificate as S
import ou3_p5_full_h_prefix_cells_v3 as H3


class Ou3P4CandidateFullWordTests(unittest.TestCase):
    def setUp(self):
        self.domain_path = C.DEFAULT_DOMAIN.resolve()
        self.domain = json.loads(self.domain_path.read_text(encoding="utf-8"))
        H3._install_backend()

    def test_30deg_ball_cover_avoids_single_cube_sqrt3_inflation(self):
        entrance = E.build(self.domain_path)
        q = float(entrance["P4_complete_word_search"]["candidate_rows"][0]["cayley_norm_upper"])
        outer = float(S.build(self.domain_path)["design_cayley_norm_upper"])
        cover = C._ball_box_cover(q, 1.5)
        self.assertEqual(len(cover), 32)
        self.assertTrue(all(C._norm_bounds_box(b)[0] <= q for b in cover))
        cover_q = max(C._norm_bounds_box(b)[1] for b in cover)
        self.assertLess(cover_q, outer)
        self.assertLess(cover_q, math.sqrt(3.0) * q)

    def test_position_entry_consumes_half_Hs_not_legacy_20m_box(self):
        C._configure_mode("H")
        e, ba_cap, meta = C._initial_error("H", self.domain)
        self.assertIsNone(ba_cap)
        self.assertAlmostEqual(meta["Hs_upper_m"], 8.5)
        self.assertAlmostEqual(meta["position_component_abs_upper_m"], 4.25)
        self.assertEqual(meta["legacy_P1_position_norm_upper_m_not_used_as_P4_entry"], 20.0)
        for i in C.H.P:
            self.assertLessEqual(e[i].abs_upper(), math.nextafter(4.25, math.inf))

    def test_A_mode_is_21_state_shipping_bias_block(self):
        C._configure_mode("A")
        src = C.H._source_cell()
        F, Q, _, ba = C._transition_and_Q("A", src, self.domain)
        self.assertEqual(C.H.N, 21)
        self.assertEqual(list(C.H.BA), [18, 19, 20])
        self.assertIsNotNone(ba)
        self.assertTrue(0.0 < ba["phi_interval"][0] <= ba["phi_interval"][1] <= 1.0)
        self.assertGreater(ba["Qd_variance_interval"][0], 0.0)
        for i in C.H.BA:
            self.assertGreater(F[i][i].lo, 0.0)
            self.assertGreater(Q[i][i].lo, 0.0)
        Hacc = C._H_acc("A", self.domain)
        for ax, i in enumerate(C.H.BA):
            self.assertEqual(Hacc[ax][i].lo, 1.0)
            self.assertEqual(Hacc[ax][i].hi, 1.0)

    def test_zero_sample_H_and_A_preflight_are_inside_outer_sector(self):
        entrance = E.build(self.domain_path)
        row = entrance["P4_complete_word_search"]["candidate_rows"][0]
        q = float(row["cayley_norm_upper"])
        outer = float(S.build(self.domain_path)["design_cayley_norm_upper"])
        box = C._ball_box_cover(q, 1.5)[0]
        for mode, dim in (("H", 18), ("A", 21)):
            r = C._run_cell(mode, self.domain_path, self.domain, box, q, outer, 0)
            self.assertEqual(r["dimension"], dim)
            self.assertTrue(r["prefix_safe_in_outer_sector"])
            self.assertLess(r["max_prefix_q_upper"], outer)

    def test_build_never_promotes_prefix_safety_to_P4_dissipation(self):
        d = C.build(self.domain_path, max_samples=0, candidate_index=0)
        self.assertEqual(C.validate(d), [])
        self.assertFalse(d["old_q8_chart_used"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])
        self.assertFalse(d["P4_NORMALIZED_CROSS_BLOCK_ESTABLISHED_HERE"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertEqual(d["H_dimension"], 18)
        self.assertEqual(d["A_dimension"], 21)


if __name__ == "__main__":
    unittest.main()
