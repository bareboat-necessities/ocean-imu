#!/usr/bin/env python3
"""Copyright 2026, Mikhail Grushinskiy. Guard external calibration temperature plumbing."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
checked = []
for sketch in sorted((root / "sensors").rglob("*.ino")):
    text = sketch.read_text()
    if not re.search(r"apply(?:Accel|Gyro)\(", text):
        continue
    assert not re.search(r"(?:std::)?isfinite\(s\.tempC\)\s*\?", text), sketch
    # These production paths share the unmodified measurement validity.
    assert "const float tempC = s.tempC;" in text, sketch
    assert re.search(r"applyAccel\([^;]+,\s*tempC\)", text), sketch
    assert re.search(r"applyGyro\([^;]+,\s*tempC\)", text), sketch
    checked.append(sketch.relative_to(root))
assert len(checked) == 6, f"Review new external calibration consumers: {checked}"
compass = (root / "src/AtomS3R/AtomS3R_CompassAppBase.h").read_text()
assert "runtime_.applyAccel(s.a, s.tempC)" in compass
assert "runtime_.applyGyro(s.w, s.tempC)" in compass
# Fusion's intentional reference temperature is NOT the sensor temperature.
# Removing these would double-apply TFG's internal thermal model.
tfg = (root / "sensors/full_marine_ins/atomS3R_ins_tfg/atomS3R_ins_tfg.ino").read_text()
assert "fusion_.update(dt_, w_cal_, a_cal_, 35.0f);" in tfg
assert "updateWaveDirection_(q_bw, 35.0f, dt_);" in tfg
print(f"test_sketch_temperature: {len(checked)} INS/diagnostic sketches + compass preserve temperature validity; TFG reference arguments intact")
