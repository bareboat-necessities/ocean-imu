#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_covariance_upper as U
import ou3_source_reachable_matrix_p3 as BASE


class CorrelatedTranslationCovarianceUpperTests(unittest.TestCase):
    def test_representative_same_history_certificate(self):
        d = U.build()
        self.assertEqual(U.validate(d), [])
        self.assertTrue(d["P2_correlation_interface_consumed"])
        self.assertEqual(d["P2_correlation_interface_version"], CORR.INTERFACE_VERSION)
        self.assertTrue(d["same_history_sufficient_statistics_used"])
        self.assertFalse(d["independent_cartesian_tau_sigma_R_S_extrema_used"])
        self.assertTrue(d["endpoint_referenced_translation_observability_used"])
        self.assertFalse(d["post_reconstruction_forward_propagation_used"])
        self.assertTrue(d["full_word_process_noise_dominator_retained"])
        self.assertFalse(d["P3_PROMOTED"])

    def test_summary_uses_one_legal_history(self):
        rt = CORR.runtime()
        start, trans = U._representative_history(rt, 137, 3, 21)
        segs = U._path_segments(start, trans, rt)
        summary = U.summarize_segments(segs, BASE.source_schedule())
        self.assertTrue(summary["all_statistics_from_one_legal_P2_history"])
        self.assertFalse(summary["independent_global_source_extrema_used"])
        self.assertEqual(summary["segments"], 3)
        self.assertGreater(summary["q_c_upper"], 0.0)
        self.assertGreater(summary["sigma_squared_upper"], 0.0)
        self.assertGreater(summary["S_measurement_variance_upper"], 0.0)

    def test_gap_theorem_fails_closed_without_progress_preserving_retarget(self):
        rt = CORR.runtime()
        start, trans = U._representative_history(rt, 137, 3, 21)
        summary = U.summarize_segments(U._path_segments(start, trans, rt), BASE.source_schedule())
        sched = dict(BASE.source_schedule())
        sched["pseudo_period_retarget_progress_preserving"] = False
        with self.assertRaisesRegex(RuntimeError, "progress-preserving"):
            U.translation_upper_from_summary(summary, 1.0, sched, require_history_cover=False)

    def test_endpoint_reference_does_not_forward_propagate_reconstruction(self):
        rt = CORR.runtime()
        start, trans = U._representative_history(rt, 729, 25, 26)
        summary = U.summarize_segments(U._path_segments(start, trans, rt), BASE.source_schedule())
        upper, timing = U.translation_upper_from_summary(
            summary, 1.0, BASE.source_schedule(), require_history_cover=True
        )
        self.assertEqual(len(upper), 4)
        self.assertTrue(all(x > 0.0 for x in upper))
        self.assertEqual(timing["translation_reference"], "word_endpoint")
        self.assertTrue(timing["endpoint_referenced_observability"])
        self.assertTrue(timing["endpoint_p_sign_similarity_applied"])
        self.assertFalse(timing["forward_propagation_after_endpoint_reconstruction"])
        self.assertTrue(timing["full_word_process_noise_dominator_retained"])
        self.assertGreater(timing["word_horizon_s_upper"], timing["observation_span_s_upper"])

    def test_constant_history_reduces_to_same_monotone_source_quantities(self):
        rt = CORR.runtime()
        # Frozen-clock-like comparison uses repeated same-cell segment kernels;
        # the sufficient statistics must not invent extrema outside that cell.
        s = 137
        segs = [CORR.segment_kernel(s, 21, rt) for _ in range(3)]
        summary = U.summarize_segments(segs, BASE.source_schedule())
        node = rt["nodes"][s]
        self.assertAlmostEqual(summary["sigma_squared_upper"], node["sigma_filter_committed_mps2"][1] ** 2, places=12)
        self.assertEqual(summary["source_nodes"], [s, s, s])


if __name__ == "__main__":
    unittest.main()