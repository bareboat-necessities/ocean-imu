#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p2_correlation_path_memory as P2
import ou3_p3_tau_decay_budget as DECAY


class P2CorrelationPathMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = P2.build()
        cls.rt = P2.runtime()

    def test_interface_is_source_only_and_frozen(self):
        d = self.payload
        self.assertEqual(P2.validate(d), [])
        self.assertEqual(d["interface_version"], "OU3_P2_CORRELATED_STAGE_TRANSFER_V1")
        self.assertEqual(d["P2_CORRELATION_INTERFACE_CERTIFICATE"], "PASS")
        self.assertTrue(d["P2_TIMING_SOURCE_MATHEMATICS_RETAINED"])
        self.assertTrue(d["tau_sigma_R_S_joint_cell_retained_per_segment"])
        self.assertTrue(d["EMA_stage_commit_history_retained"])
        self.assertFalse(d["arbitrary_cartesian_tuner_switching_used"])
        self.assertFalse(d["old_800_node_ancestor_hull_allowed_for_correlated_P3"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["filter_changed"])

    def test_pair_state_and_gap_semantics_are_exact(self):
        d = self.payload
        self.assertEqual(d["physical_source_states"], 800)
        self.assertEqual(d["clock_gap_alphabet_samples"], list(range(13, 27)))
        self.assertGreater(d["stage_boundary_pair_states"], 0)
        self.assertGreater(d["gap_labelled_first_order_edges"], 0)
        self.assertGreater(d["gap_labelled_pair_state_edges_factorized_count"], 0)

        # Pick the first legal pair and one legal labelled continuation.  The
        # transition must shift (c,s) -> (s,t) without losing the staged cell.
        c = next(i for i, out in enumerate(self.rt["union_successors"]) if out)
        s = min(self.rt["union_successors"][c])
        gap = next(g for g in self.rt["gaps"] if P2.successors(s, g, self.rt))
        t = P2.successors(s, gap, self.rt)[0]
        tr = P2.transition(c, s, gap, t, self.rt)
        self.assertEqual(tr["start_pair"], [c, s])
        self.assertEqual(tr["end_pair"], [s, t])
        self.assertEqual(tr["boundary_sample_applied_node"], c)
        self.assertEqual(tr["following_segment_applied_node"], s)
        self.assertTrue(tr["pair_shift_exact"])
        self.assertTrue(tr["same_staged_node_becomes_next_committed"])

    def test_segment_statistics_keep_one_common_source_node(self):
        for s in (0, 137, 729, 799):
            for gap in (13, 21, 26):
                k = P2.segment_kernel(s, gap, self.rt)
                node = self.rt["nodes"][s]
                self.assertEqual(k["applied_source_node"], s)
                self.assertEqual(k["samples"], gap)
                self.assertEqual(k["tau_s"], node["tau_s"])
                self.assertEqual(k["sigma_filter_committed_mps2"], node["sigma_filter_committed_mps2"])
                self.assertEqual(k["R_S_filter_std"], node["R_S_filter_std"])
                self.assertTrue(k["tau_sigma_R_S_from_same_physical_cell"])
                self.assertFalse(k["independent_coordinate_extrema_used"])
                for key in (
                    "duration_s", "lambda_inv_tau_per_s", "decay_exponent_integral",
                    "q_c_m2ps5", "q_c_time_mass", "sigma_squared",
                    "inverse_R_S_variance", "normalized_S_information_per_packet",
                ):
                    lo, hi = map(float, k[key])
                    self.assertTrue(math.isfinite(lo) and math.isfinite(hi))
                    self.assertGreater(lo, 0.0)
                    self.assertGreaterEqual(hi, lo)

    def test_consumer_contract_forbids_flat_correlated_extrema(self):
        c = self.payload["consumer_contract"]
        self.assertTrue(c["correlated_quantities_must_come_from_same_segment_node"])
        self.assertTrue(c["legal_word_must_follow_pair_shift_transition_rule"])
        self.assertEqual(c["independent_tau_sigma_R_S_extremization_before_propagation"], "FORBIDDEN")
        self.assertEqual(c["global_800_ancestor_hull_as_P3_covariance_information_input"], "FORBIDDEN")

    def test_tau_decay_is_explicit_permitted_scalar_projection_consumer(self):
        d = DECAY.build(window_samples=(50,))
        self.assertEqual(DECAY.validate(d), [])
        self.assertEqual(d["P2_correlation_interface_version"], P2.INTERFACE_VERSION)
        self.assertTrue(d["P2_correlation_interface_consumed"])
        self.assertTrue(d["projection_permitted_by_P2_contract"])
        self.assertEqual(d["projection_role"], "TAU_ONLY_SCALAR_MAXIMUM")
        self.assertTrue(d["sigma_RS_projection_only_adds_paths"])
        self.assertFalse(d["flat_800_node_ancestor_hull_consumed"])

    def test_tau_consumer_does_not_rebuild_old_flat_source_graph(self):
        text = (TOOLS / "ou3_p3_tau_decay_budget.py").read_text(encoding="utf-8")
        self.assertIn("import ou3_p2_correlation_path_memory as CORR", text)
        self.assertNotIn("import ou3_p4_sample_clock_source_refinement", text)
        self.assertNotIn("import ou3_p4_source_path_reachability", text)


if __name__ == "__main__":
    unittest.main()
