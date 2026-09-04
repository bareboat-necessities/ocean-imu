#!/usr/bin/env python3
"""Temporary branch helper for the scheduler-progress proof cleanup."""
from pathlib import Path

# Migrate the ordered point diagnostic from the historical fmod setter to the
# deployed progress-preserving retarget and consume the positive recurrence
# certificate.
p = Path("tools/ou3_p3_ordered_witness_covariance_diagnostic.py")
text = p.read_text(encoding="utf-8")
repls = [
    (
        "import ou3_p3_pseudo_scheduler_starvation_witness as TIMER\n",
        "import ou3_p3_pseudo_scheduler_progress_certificate as TIMER\n",
    ),
    (
        '''def _set_period(elapsed, period):\n    """Use the exact binary32 setter transcription shared with the scheduler witness."""\n    return TIMER._set_period(elapsed, period)\n''',
        '''def _set_period(elapsed, period):\n    """Use the exact binary32 progress-preserving shipping retarget."""\n    return TIMER._retarget(elapsed, period)\n''',
    ),
    (
        '''    rt = CORR.runtime(path); sched = BASE.source_schedule(); h = TIMER._f32(rt["clock"]["dt_binary32_s"])\n    target = HIST._global_word_target(domain, sched, h); N = int(target["target_samples"])\n''',
        '''    rt = CORR.runtime(path); sched = BASE.source_schedule(); h = TIMER._f32(rt["clock"]["dt_binary32_s"])\n    progress = TIMER.build(path)\n    progress_failures = TIMER.validate(progress)\n    if progress_failures:\n        raise RuntimeError(f"scheduler progress prerequisite failed: {progress_failures}")\n    target = HIST._global_word_target(domain, sched, h); N = int(target["target_samples"])\n''',
    ),
    (
        '''        "exact_witness_source_order_retained": True, "exact_gap_labelled_legal_extension_used": True,\n        "pseudo_period_change_uses_fmod_semantics": True, "periodic_update_due_shipping_semantics_transcribed": True,\n        "pseudo_scheduler_numeric_type": "binary32/float",\n''',
        '''        "exact_witness_source_order_retained": True, "exact_gap_labelled_legal_extension_used": True,\n        "pseudo_period_change_preserves_elapsed_service_credit": True, "pseudo_period_change_uses_fmod_semantics": False,\n        "periodic_update_due_shipping_semantics_transcribed": True,\n        "scheduler_progress_certificate_consumed": True,\n        "scheduler_uniform_max_gap_samples": int(progress["certified_uniform_max_gap_samples"]),\n        "pseudo_scheduler_numeric_type": "binary32/float",\n''',
    ),
    (
        '''        "next_obligation": "select a certified time-ordered covariance-upper quotient carrying source order and pseudo-scheduler phase; do not promote from this point diagnostic",\n''',
        '''        "next_obligation": "rerun the canonical source-reachable P3 covariance closure under the certified progress-preserving S recurrence; do not promote from this point diagnostic",\n''',
    ),
    (
        '''    for k in ("diagnostic_only", "P2_correlation_interface_consumed", "exact_witness_source_order_retained",\n              "exact_gap_labelled_legal_extension_used", "pseudo_period_change_uses_fmod_semantics",\n              "periodic_update_due_shipping_semantics_transcribed", "source_cells_use_one_real_upper_corner"):\n''',
        '''    for k in ("diagnostic_only", "P2_correlation_interface_consumed", "exact_witness_source_order_retained",\n              "exact_gap_labelled_legal_extension_used", "pseudo_period_change_preserves_elapsed_service_credit",\n              "periodic_update_due_shipping_semantics_transcribed", "scheduler_progress_certificate_consumed",\n              "source_cells_use_one_real_upper_corner"):\n''',
    ),
    (
        '''    for k in ("trajectory_replay_used", "filter_changed", "declared_domain_changed", "canonical_gate_changed",\n              "accelerometer_measurement_updates_credited", "interval_certificate", "uniform_covariance_upper_certificate",\n              "matched_margin_computed", "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED"):\n''',
        '''    for k in ("trajectory_replay_used", "filter_changed", "declared_domain_changed", "canonical_gate_changed",\n              "pseudo_period_change_uses_fmod_semantics", "accelerometer_measurement_updates_credited",\n              "interval_certificate", "uniform_covariance_upper_certificate", "matched_margin_computed",\n              "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED"):\n''',
    ),
    (
        '''    if d.get("pseudo_scheduler_numeric_type") != "binary32/float": f.append("ordered diagnostic is not using shipping float scheduler arithmetic")\n''',
        '''    if d.get("pseudo_scheduler_numeric_type") != "binary32/float": f.append("ordered diagnostic is not using shipping float scheduler arithmetic")\n    if int(d.get("scheduler_uniform_max_gap_samples", 0)) != 33: f.append("ordered diagnostic lost the 33-sample scheduler recurrence")\n''',
    ),
]
for old, new in repls:
    if text.count(old) != 1:
        raise SystemExit(f"ordered diagnostic replacement anchor count={text.count(old)} for {old[:80]!r}")
    text = text.replace(old, new, 1)
p.write_text(text, encoding="utf-8")

# Replace the ordered-diagnostic unit test that asserted the retired fmod setter.
Path("tests/validation/test_ou3_p3_ordered_witness_covariance_diagnostic.py").write_text('''from pathlib import Path\nimport sys\nimport unittest\n\nsys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))\n\nimport ou3_p3_ordered_witness_covariance_diagnostic as D\n\n\nclass OrderedWitnessCovarianceDiagnosticTests(unittest.TestCase):\n    def test_period_change_preserves_unexpired_service_credit(self):\n        elapsed = D.TIMER._f32(0.12)\n        period = D.TIMER._f32(0.14)\n        self.assertEqual(D._set_period(elapsed, period), elapsed)\n\n    def test_overdue_period_change_arms_service(self):\n        period = D.TIMER._f32(0.13)\n        armed = D._set_period(D.TIMER._f32(0.14), period)\n        self.assertEqual(armed, D.TIMER._nextafterf_down_positive(period))\n        due, _ = D._due(D.TIMER._f32(0.005), period, armed)\n        self.assertTrue(due)\n\n    def test_periodic_due_waits_then_fires(self):\n        due, elapsed = D._due(0.005, 0.015, 0.0)\n        self.assertFalse(due)\n        due, elapsed = D._due(0.005, 0.015, elapsed)\n        self.assertFalse(due)\n        due, elapsed = D._due(0.005, 0.015, elapsed)\n        self.assertTrue(due)\n        self.assertGreaterEqual(elapsed, 0.0)\n        self.assertLess(elapsed, 0.015)\n\n    def test_extension_prefers_exact_self_edge_and_ends_on_target_phase(self):\n        rt = {\n            "nodes": [{}, {}],\n            "gaps": [13, 14],\n            "labelled_successors": [\n                [[0], [1]],\n                [[1], [0]],\n            ],\n        }\n        witness = [{\n            "source": 0,\n            "successor": 0,\n            "gap_samples": 13,\n            "cumulative_samples": 13,\n        }]\n        segs = D.extend_witness_to_target(witness, rt, 30)\n        self.assertEqual([row[0] for row in segs], [0, 0, 0])\n        self.assertEqual([row[1] for row in segs], [13, 13, 4])\n        self.assertEqual(segs[-1][2], 13)\n        self.assertFalse(segs[-1][4])\n\n    def test_synthetic_tuple_reconstructs_process_maxima(self):\n        summary = {\n            "sigma_squared_upper": 36.0,\n            "q_c_upper": 72.0,\n            "S_measurement_variance_upper": 16.0,\n            "pseudo_update_cadence_s": [0.005, 0.25],\n        }\n        p = D._synthetic(summary)\n        self.assertAlmostEqual(p["sigma"] ** 2, 36.0)\n        self.assertAlmostEqual(2.0 * p["sigma"] ** 2 / p["tau"], 72.0)\n        self.assertEqual(p["period"], 0.25)\n        self.assertEqual(p["Rstd"], 4.0)\n        self.assertIsNone(p["node"])\n\n\nif __name__ == "__main__":\n    unittest.main()\n''', encoding="utf-8")

# Historical zero-fire producers are retained in git history (#476) but are no
# longer live shipping diagnostics after the setter repair.
for dead in (
    "tools/ou3_p3_pseudo_scheduler_starvation_witness.py",
    "tests/validation/test_ou3_p3_pseudo_scheduler_starvation_witness.py",
    "tools/ou3_p3_tau_ema_scheduler_cycle_diagnostic.py",
    "tests/validation/test_ou3_p3_tau_ema_scheduler_cycle_diagnostic.py",
):
    q = Path(dead)
    if not q.exists():
        raise SystemExit(f"expected historical diagnostic missing: {dead}")
    q.unlink()
