#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for relative in [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]:
    path = ROOT / relative
    text = path.read_text(encoding="utf-8")
    old = "            dt, dir_filter_.getLastStableConfidence());\n"
    new = "            dt, dir_filter_.getConfidence());\n"
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{relative}: expected one confidence gate, found {count}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")

print("Applied current axis-confidence gating to OU-II and OU-III.")
