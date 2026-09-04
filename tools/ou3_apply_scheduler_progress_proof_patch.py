#!/usr/bin/env python3
"""Temporary branch helper: bind the canonical P3 schedule to scheduler progress."""
from pathlib import Path

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
