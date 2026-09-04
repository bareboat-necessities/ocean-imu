#!/usr/bin/env python3
"""Temporary branch helper: bind every retained P3 gap theorem to scheduler progress.

This helper is deliberately idempotent because an earlier staging pass may have
already committed the source-uniform binding before a later pass adds the
same-history binding.  Already-bound files are verified, not rewritten.
"""
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{label} anchor count={n}")
    return text.replace(old, new, 1)


# Canonical source schedule and source-uniform translation upper.
p = Path("tools/ou3_source_reachable_matrix_p3.py")
text = p.read_text(encoding="utf-8")
if '"pseudo_period_retarget_progress_preserving":True' not in text:
    text = replace_once(
        text,
        '''REPO = Path(__file__).resolve().parents[1]\nWRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"\nDEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"\n''',
        '''REPO = Path(__file__).resolve().parents[1]\nWRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"\nMEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"\nCORE_MATH = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"\nDEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"\n''',
        "path",
    )
    text = replace_once(
        text,
        '''def source_schedule() -> dict:\n    text = WRAPPER.read_text(encoding="utf-8")\n    for marker in (\n''',
        '''def source_schedule() -> dict:\n    text = WRAPPER.read_text(encoding="utf-8")\n    mekf_text = MEKF.read_text(encoding="utf-8")\n    core_text = CORE_MATH.read_text(encoding="utf-8")\n    for marker in (\n''',
        "source_schedule",
    )
    text = replace_once(
        text,
        '''        if marker not in text:\n            raise RuntimeError(f"missing deployed schedule semantic: {marker}")\n    names = (\n''',
        '''        if marker not in text:\n            raise RuntimeError(f"missing deployed schedule semantic: {marker}")\n    for marker in (\n        "retarget_period_elapsed_progress_preserving(",\n        "pseudo_update_elapsed_s_, new_period",\n        "pseudo_update_period_s_ = new_period;",\n    ):\n        if marker not in mekf_text:\n            raise RuntimeError(f"missing progress-preserving OU-III scheduler semantic: {marker}")\n    for marker in (\n        "if (elapsed < period) return elapsed;",\n        "return std::nextafter(period, T(0));",\n    ):\n        if marker not in core_text:\n            raise RuntimeError(f"missing progress-preserving scheduler helper semantic: {marker}")\n    if "std::fmod(pseudo_update_elapsed_s_, pseudo_update_period_s_)" in mekf_text:\n        raise RuntimeError("legacy pseudo-period setter still discards elapsed service credit")\n    names = (\n''',
        "schedule markers",
    )
    text = replace_once(
        text,
        '''        "dt_s":c["FREQ_SMOOTHER_DT"],\n        "proof_kind":"SOURCE_REACHABLE_INVARIANT_CELL_OVERAPPROXIMATION",\n''',
        '''        "dt_s":c["FREQ_SMOOTHER_DT"],\n        "pseudo_period_retarget_progress_preserving":True,\n        "proof_kind":"SOURCE_REACHABLE_INVARIANT_CELL_OVERAPPROXIMATION",\n''',
        "schedule return",
    )
    text = replace_once(
        text,
        '''def translation_upper(tau: Interval,sigma: Interval,rs: Interval,Tpe: float,sched: dict) -> tuple[list[float],dict]:\n    h=sched["dt_s"]\n    cadence=cadence_bounds(tau,sched)\n''',
        '''def translation_upper(tau: Interval,sigma: Interval,rs: Interval,Tpe: float,sched: dict) -> tuple[list[float],dict]:\n    h=sched["dt_s"]\n    if sched.get("pseudo_period_retarget_progress_preserving") is not True:\n        raise RuntimeError("finite S-observation gap requires progress-preserving pseudo-period retargeting")\n    cadence=cadence_bounds(tau,sched)\n''',
        "translation upper",
    )
    p.write_text(text, encoding="utf-8")
else:
    for marker in (
        'MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"',
        'CORE_MATH = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"',
        '"retarget_period_elapsed_progress_preserving("',
        '"pseudo_period_retarget_progress_preserving":True',
        'finite S-observation gap requires progress-preserving pseudo-period retargeting',
    ):
        if marker not in text:
            raise SystemExit(f"already-bound source backend lost marker: {marker}")


# Same-history upper used by the actual canonical P2-V1 history frontier.
p = Path("tools/ou3_p3_correlated_translation_covariance_upper.py")
text = p.read_text(encoding="utf-8")
if '"progress_preserving_scheduler_required_for_gap_bound": True' not in text:
    text = replace_once(
        text,
        '''* every S-observation gap is bounded by the path maximum cadence plus one sample;\n''',
        '''* with progress-preserving pseudo-period retargeting, every S-observation gap is bounded by the path maximum cadence plus one sample;\n''',
        "same-history doc",
    )
    text = replace_once(
        text,
        '''    if summary.get("independent_global_source_extrema_used") is not False:\n        raise ValueError("independent global source extrema are forbidden")\n\n    h = float(sched["dt_s"])\n''',
        '''    if summary.get("independent_global_source_extrema_used") is not False:\n        raise ValueError("independent global source extrema are forbidden")\n    if sched.get("pseudo_period_retarget_progress_preserving") is not True:\n        raise RuntimeError(\n            "same-history finite S-observation gap requires progress-preserving pseudo-period retargeting"\n        )\n\n    h = float(sched["dt_s"])\n''',
        "same-history schedule",
    )
    text = replace_once(
        text,
        '''        "retained_translation_observability_theorem_reused": True,\n        "monotone_path_maxima_only": True,\n''',
        '''        "retained_translation_observability_theorem_reused": True,\n        "progress_preserving_scheduler_required_for_gap_bound": True,\n        "monotone_path_maxima_only": True,\n''',
        "same-history build flag",
    )
    text = replace_once(
        text,
        '''        "same_history_sufficient_statistics_used",\n        "retained_translation_observability_theorem_reused",\n        "monotone_path_maxima_only", "full_covariance_word_history_required",\n''',
        '''        "same_history_sufficient_statistics_used",\n        "retained_translation_observability_theorem_reused",\n        "progress_preserving_scheduler_required_for_gap_bound",\n        "monotone_path_maxima_only", "full_covariance_word_history_required",\n''',
        "same-history validate flag",
    )
    p.write_text(text, encoding="utf-8")
else:
    for marker in (
        "same-history finite S-observation gap requires progress-preserving pseudo-period retargeting",
        '"progress_preserving_scheduler_required_for_gap_bound": True',
    ):
        if marker not in text:
            raise SystemExit(f"already-bound same-history backend lost marker: {marker}")


# Regression: removing the schedule contract must make the same-history theorem
# refuse to produce a finite observation-gap upper.
p = Path("tests/validation/test_ou3_p3_correlated_translation_covariance_upper.py")
text = p.read_text(encoding="utf-8")
method = "def test_gap_theorem_fails_closed_without_progress_preserving_retarget(self):"
if method not in text:
    anchor = '''    def test_constant_history_reduces_to_same_monotone_source_quantities(self):\n'''
    block = '''    def test_gap_theorem_fails_closed_without_progress_preserving_retarget(self):\n        rt = CORR.runtime()\n        start, trans = U._representative_history(rt, 137, 3, 21)\n        summary = U.summarize_segments(U._path_segments(start, trans, rt), BASE.source_schedule())\n        sched = dict(BASE.source_schedule())\n        sched["pseudo_period_retarget_progress_preserving"] = False\n        with self.assertRaisesRegex(RuntimeError, "progress-preserving"):\n            U.translation_upper_from_summary(summary, 1.0, sched, require_history_cover=False)\n\n    def test_constant_history_reduces_to_same_monotone_source_quantities(self):\n'''
    text = replace_once(text, anchor, block, "same-history test")
    p.write_text(text, encoding="utf-8")
