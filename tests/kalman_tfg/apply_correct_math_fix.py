#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "src/kalman_tfg/Kalman3D_Wave_TFG.h"
text = path.read_text()
old = "    void set_accel_bias_limit(T limit) { acc_bias_limit_ = std::abs(limit); }\n"
if old in text:
    text = text.replace(old, "", 1)
path.write_text(text)
print("Removed obsolete accelerometer-bias clamp setter")
