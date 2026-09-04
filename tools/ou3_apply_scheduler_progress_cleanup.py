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

# Replace the dedicated structural workflow with the positive scheduler
# recurrence + ordered point diagnostic.  Artifacts upload even on failure.
Path(".github/workflows/ou3-p3-whole-word-probe.yml").write_text(r'''name: ou3-p3-scheduler-progress-and-ordered-covariance

on:
  pull_request:
    paths:
      - "tools/ou3_p3_four_max_global_label_witness.py"
      - "tests/validation/test_ou3_p3_four_max_global_label_witness.py"
      - "tools/ou3_p3_pseudo_scheduler_progress_certificate.py"
      - "tests/validation/test_ou3_p3_pseudo_scheduler_progress_certificate.py"
      - "tools/ou3_p3_ordered_witness_covariance_diagnostic.py"
      - "tests/validation/test_ou3_p3_ordered_witness_covariance_diagnostic.py"
      - ".github/workflows/ou3-p3-whole-word-probe.yml"
      - "tools/ou3_interval.py"
      - "tools/ou3_p2_*.py"
      - "tools/ou3_p3_*.py"
      - "tools/ou3_p4_*.py"
      - "tools/ou3_source_reachable_matrix_p3.py"
      - "tools/ou3_proof_operating_domain.json"
      - "src/kalman_ou_iii/**"
      - "src/kalman_ou_common/**"
      - "src/tuner/**"
  workflow_dispatch:

permissions:
  contents: read

jobs:
  scheduler-progress-and-ordered-covariance:
    runs-on: ubuntu-latest
    timeout-minutes: 120
    steps:
      - uses: actions/checkout@v6
        with:
          persist-credentials: false

      - name: Compile live P3 scheduler and ordered diagnostics
        run: |
          python3 -m py_compile \
            tools/ou3_p3_four_max_global_label_witness.py \
            tools/ou3_p3_pseudo_scheduler_progress_certificate.py \
            tools/ou3_p3_ordered_witness_covariance_diagnostic.py \
            tests/validation/test_ou3_p3_four_max_global_label_witness.py \
            tests/validation/test_ou3_p3_pseudo_scheduler_progress_certificate.py \
            tests/validation/test_ou3_p3_ordered_witness_covariance_diagnostic.py

      - name: Test witness, scheduler progress, and ordered extension invariants
        working-directory: tests/validation
        run: |
          python3 -m unittest -v \
            test_ou3_p3_four_max_global_label_witness \
            test_ou3_p3_pseudo_scheduler_progress_certificate \
            test_ou3_p3_ordered_witness_covariance_diagnostic

      - name: Certify progress-preserving pseudo-scheduler recurrence
        run: |
          python3 tools/ou3_p3_pseudo_scheduler_progress_certificate.py \
            --output /tmp/ou3_p3_pseudo_scheduler_progress_certificate.json

      - name: Enforce scheduler recurrence contract
        run: |
          python3 - <<'PY'
          import json
          d=json.load(open('/tmp/ou3_p3_pseudo_scheduler_progress_certificate.json'))
          assert d['validation_pass'], d['validation_failures']
          assert d['core_progress_helper_markers_present'] is True
          assert d['ou3_setter_progress_markers_present'] is True
          assert d['legacy_setter_fmod_progress_reset_absent'] is True
          assert d['shipping_scheduler_numeric_type'] == 'binary32/float'
          assert d['fixed_max_period_first_fire_samples'] == 33
          assert d['certified_uniform_max_gap_samples'] == 33
          assert d['uniform_gap_within_translation_upper_bound'] is True
          assert d['former_635_sample_zero_fire_cycle_broken'] is True
          assert d['former_starvation_cycle_pseudo_firings'] > 0
          assert d['former_starvation_cycle_worst_gap_samples'] <= 33
          assert d['scheduler_recurrence_certificate'] is True
          assert d['P3_PROMOTED'] is False
          assert d['P4_PROMOTED'] is False
          assert d['P5_PROMOTED'] is False
          print('PERIOD_RANGE', d['deployed_pseudo_period_binary32_s'])
          print('MAX_GAP_SAMPLES', d['certified_uniform_max_gap_samples'])
          print('MAX_GAP_S', d['certified_uniform_max_gap_s'])
          print('TRANSLATION_BOUND_S', d['translation_upper_cadence_plus_h_s'])
          print('FORMER_CYCLE_FIRINGS', d['former_starvation_cycle_pseudo_firings'])
          print('FORMER_CYCLE_WORST_GAP', d['former_starvation_cycle_worst_gap_samples'])
          PY

      - name: Find shortest legal global four-max label witness
        run: |
          python3 tools/ou3_p3_four_max_global_label_witness.py \
            --output /tmp/ou3_p3_four_max_global_label_witness.json

      - name: Enforce four-max collapse witness contract
        run: |
          python3 - <<'PY'
          import json
          d=json.load(open('/tmp/ou3_p3_four_max_global_label_witness.json'))
          assert d['validation_pass'], d['validation_failures']
          assert d['target_samples'] == 635
          assert d['global_label_reachable_within_word'] is True
          assert d['P3_PROMOTED'] is False
          assert d['P4_PROMOTED'] is False
          assert d['P5_PROMOTED'] is False
          print('GLOBAL_LABEL', d['global_adverse_rank_label'])
          print('MINIMUM_COST_SAMPLES', d['global_label_minimum_cost_samples'])
          PY

      - name: Run ordered source and progress-preserving scheduler point diagnostic
        run: |
          python3 tools/ou3_p3_ordered_witness_covariance_diagnostic.py \
            --output /tmp/ou3_p3_ordered_witness_covariance_diagnostic.json

      - name: Enforce ordered diagnostic non-promotion contract
        run: |
          python3 - <<'PY'
          import json
          d=json.load(open('/tmp/ou3_p3_ordered_witness_covariance_diagnostic.json'))
          assert d['validation_pass'], d['validation_failures']
          assert d['target_samples'] == 635
          assert d['exact_witness_source_order_retained'] is True
          assert d['exact_gap_labelled_legal_extension_used'] is True
          assert d['pseudo_period_change_preserves_elapsed_service_credit'] is True
          assert d['pseudo_period_change_uses_fmod_semantics'] is False
          assert d['scheduler_progress_certificate_consumed'] is True
          assert d['scheduler_uniform_max_gap_samples'] == 33
          assert d['periodic_update_due_shipping_semantics_transcribed'] is True
          assert d['accelerometer_measurement_updates_credited'] is False
          assert d['interval_certificate'] is False
          assert d['uniform_covariance_upper_certificate'] is False
          assert d['P3_PROMOTED'] is False
          assert d['P4_PROMOTED'] is False
          assert d['P5_PROMOTED'] is False
          print('OLD_UPPER_STD', d['old_four_max_upper_std'])
          print('ORDERED_STD', d['ordered_phase_envelope_std'])
          print('ORDERED_VARIANCE_GAIN', d['old_upper_to_ordered_variance_gain'])
          print('SYNTHETIC_STD', d['synthetic_phase_envelope_std'])
          print('SYNTHETIC_VARIANCE_GAIN', d['old_upper_to_synthetic_variance_gain'])
          print('NEXT', d['next_obligation'])
          PY

      - name: Upload scheduler progress certificate artifact
        if: always()
        uses: actions/upload-artifact@v6
        with:
          name: ou3-p3-pseudo-scheduler-progress-certificate
          path: /tmp/ou3_p3_pseudo_scheduler_progress_certificate.json
          retention-days: 7

      - name: Upload four-max witness artifact
        if: always()
        uses: actions/upload-artifact@v6
        with:
          name: ou3-p3-four-max-global-label-witness
          path: /tmp/ou3_p3_four_max_global_label_witness.json
          retention-days: 7

      - name: Upload ordered witness covariance diagnostic artifact
        if: always()
        uses: actions/upload-artifact@v6
        with:
          name: ou3-p3-ordered-witness-covariance-diagnostic
          path: /tmp/ou3_p3_ordered_witness_covariance_diagnostic.json
          retention-days: 7
''', encoding="utf-8")
