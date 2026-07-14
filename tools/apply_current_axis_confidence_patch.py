#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for relative in [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]:
    path = ROOT / relative
    text = path.read_text(encoding="utf-8")

    old_frame = (
        "        const auto direction_accel = wave_direction::heading_frame_acceleration<float>(\n"
        "            mekf_->quaternion_boat(), acc, g_std);\n\n"
    )
    new_frame = old_frame + (
        "        if (!direction_accel.heading_valid) {\n"
        "            dir_sign_state_ = UNCERTAIN;\n"
        "            return;\n"
        "        }\n\n"
    )
    frame_count = text.count(old_frame)
    if frame_count != 1:
        raise RuntimeError(
            f"{relative}: expected one direction frame conversion, found {frame_count}")
    text = text.replace(old_frame, new_frame, 1)

    old_confidence = "            dt, dir_filter_.getLastStableConfidence());\n"
    new_confidence = "            dt, dir_filter_.getConfidence());\n"
    confidence_count = text.count(old_confidence)
    if confidence_count != 1:
        raise RuntimeError(
            f"{relative}: expected one confidence gate, found {confidence_count}")
    text = text.replace(old_confidence, new_confidence, 1)

    path.write_text(text, encoding="utf-8")

print("Applied valid-heading and current axis-confidence gating to OU-II and OU-III.")
