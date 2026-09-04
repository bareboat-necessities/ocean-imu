#!/usr/bin/env python3
"""Temporary branch helper: bind every retained P3 gap theorem to scheduler progress."""
from pathlib import Path

# Canonical source schedule and source-uniform translation upper.
p = Path("tools/ou3_source_reachable_matrix_p3.py")
text = p.read_text(encoding="utf-8")

old = '''REPO = Path(__file__).resolve().parents[1]\nWRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"\nDEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"\n'''
new = '''REPO = Path(__file__).resolve().parents[1]\nWRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"\nMEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"\nCORE_MATH = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"\nDEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"\n'''
if text.count(old) != 1:
    raise SystemExit(f"path anchor count={text.count(old)}")
text = text.replace(old, new, 1)

old = '''def source_schedule() -> dict:\n    text = WRAPPER.read_text(encoding="utf-8")\n    for marker in (\n'''
new = '''def source_schedule() -> dict:\n    text = WRAPPER.read_text(encoding="utf-8")\n    mekf_text = MEKF.read_text(encoding="utf-8")\n    core_text = CORE_MATH.read_text(encoding="utf-8")\n    for marker in (\n'''
if text.count(old) != 1:
    raise SystemExit(f"source_schedule anchor count={text.count(old)}")
text = text.replace(old, new, 1)

old = '''        if marker not in text:\n            raise RuntimeError(f"missing deployed schedule semantic: {marker}")\n    names = (\n'''
new = '''        if marker not in text:\n            raise RuntimeError(f"missing deployed schedule semantic: {marker}")\n    for marker in (\n        "retarget_period_elapsed_progress_preserving(",\n        "pseudo_update_elapsed_s_, new_period",\n        "pseudo_update_period_s_ = new_period;",\n    ):\n        if marker not in mekf_text:\n            raise RuntimeError(f"missing progress-preserving OU-III scheduler semantic: {marker}")\n    for marker in (\n        "if (elapsed < period) return elapsed;",\n        "return std::nextafter(period, T(0));",\n    ):\n        if marker not in core_text:\n            raise RuntimeError(f"missing progress-preserving scheduler helper semantic: {marker}")\n    if "std::fmod(pseudo_update_elapsed_s_, pseudo_update_period_s_)" in mekf_text:\n        raise RuntimeError("legacy pseudo-period setter still discards elapsed service credit")\n    names = (\n'''
if text.count(old) != 1:
    raise SystemExit(f"schedule marker anchor count={text.count(old)}")
text = text.replace(old, new, 1)

old = '''        "dt_s":c["FREQ_SMOOTHER_DT"],\n        "proof_kind":"SOURCE_REACHABLE_INVARIANT_CELL_OVERAPPROXIMATION",\n'''
new = '''        "dt_s":c["FREQ_SMOOTHER_DT"],\n        "pseudo_period_retarget_progress_preserving":True,\n        "proof_kind":"SOURCE_REACHABLE_INVARIANT_CELL_OVERAPPROXIMATION",\n'''
if text.count(old) != 1:
    raise SystemExit(f"schedule return anchor count={text.count(old)}")
text = text.replace(old, new, 1)

old = '''def translation_upper(tau: Interval,sigma: Interval,rs: Interval,Tpe: float,sched: dict) -> tuple[list[float],dict]:\n    h=sched["dt_s"]\n    cadence=cadence_bounds(tau,sched)\n'''
new = '''def translation_upper(tau: Interval,sigma: Interval,rs: Interval,Tpe: float,sched: dict) -> tuple[list[float],dict]:\n    h=sched["dt_s"]\n    if sched.get("pseudo_period_retarget_progress_preserving") is not True:\n        raise RuntimeError("finite S-observation gap requires progress-preserving pseudo-period retargeting")\n    cadence=cadence_bounds(tau,sched)\n'''
if text.count(old) != 1:
    raise SystemExit(f"translation upper anchor count={text.count(old)}")
text = text.replace(old, new, 1)
p.write_text(text, encoding="utf-8")

# Same-history upper used by the actual canonical P2-V1 history frontier. It
# must fail closed on the same implementation contract; otherwise canonical P3
# could keep using cadence_max+h after a future setter regression.
p = Path("tools/ou3_p3_correlated_translation_covariance_upper.py")
text = p.read_text(encoding="utf-8")
old = '''* every S-observation gap is bounded by the path maximum cadence plus one sample;\n'''
new = '''* with progress-preserving pseudo-period retargeting, every S-observation gap is bounded by the path maximum cadence plus one sample;\n'''
if text.count(old) != 1:
    raise SystemExit(f"same-history doc anchor count={text.count(old)}")
text = text.replace(old, new, 1)
old = '''    if summary.get("independent_global_source_extrema_used") is not False:\n        raise ValueError("independent global source extrema are forbidden")\n\n    h = float(sched["dt_s"])\n'''
new = '''    if summary.get("independent_global_source_extrema_used") is not False:\n        raise ValueError("independent global source extrema are forbidden")\n    if sched.get("pseudo_period_retarget_progress_preserving") is not True:\n        raise RuntimeError(\n            "same-history finite S-observation gap requires progress-preserving pseudo-period retargeting"\n        )\n\n    h = float(sched["dt_s"])\n'''
if text.count(old) != 1:
    raise SystemExit(f"same-history upper schedule anchor count={text.count(old)}")
text = text.replace(old, new, 1)
old = '''        "retained_translation_observability_theorem_reused": True,\n        "monotone_path_maxima_only": True,\n'''
new = '''        "retained_translation_observability_theorem_reused": True,\n        "progress_preserving_scheduler_required_for_gap_bound": True,\n        "monotone_path_maxima_only": True,\n'''
if text.count(old) != 1:
    raise SystemExit(f"same-history build flag anchor count={text.count(old)}")
text = text.replace(old, new, 1)
old = '''        "same_history_sufficient_statistics_used",\n        "retained_translation_observability_theorem_reused",\n        "monotone_path_maxima_only", "full_covariance_word_history_required",\n'''
new = '''        "same_history_sufficient_statistics_used",\n        "retained_translation_observability_theorem_reused",\n        "progress_preserving_scheduler_required_for_gap_bound",\n        "monotone_path_maxima_only", "full_covariance_word_history_required",\n'''
if text.count(old) != 1:
    raise SystemExit(f"same-history validate flag anchor count={text.count(old)}")
text = text.replace(old, new, 1)
p.write_text(text, encoding="utf-8")

# Regression: deleting the schedule contract must make the same-history theorem
# refuse to produce a finite observation-gap upper.
p = Path("tests/validation/test_ou3_p3_correlated_translation_covariance_upper.py")
text = p.read_text(encoding="utf-8")
anchor = '''    def test_constant_history_reduces_to_same_monotone_source_quantities(self):\n'''
block = '''    def test_gap_theorem_fails_closed_without_progress_preserving_retarget(self):\n        rt = CORR.runtime()\n        start, trans = U._representative_history(rt, 137, 3, 21)\n        summary = U.summarize_segments(U._path_segments(start, trans, rt), BASE.source_schedule())\n        sched = dict(BASE.source_schedule())\n        sched["pseudo_period_retarget_progress_preserving"] = False\n        with self.assertRaisesRegex(RuntimeError, "progress-preserving"):\n            U.translation_upper_from_summary(summary, 1.0, sched, require_history_cover=False)\n\n    def test_constant_history_reduces_to_same_monotone_source_quantities(self):\n'''
if text.count(anchor) != 1:
    raise SystemExit(f"same-history test anchor count={text.count(anchor)}")
text = text.replace(anchor, block, 1)
p.write_text(text, encoding="utf-8")
