from pathlib import Path

p = Path('src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h')
text = p.read_text()
old = '    void enableTuner(bool flag = true) { enable_tuner_ = flag; }\n'
new = '    void enableTuner(bool flag = true) {\n        enable_tuner_ = flag;\n    }\n'
if old not in text:
    raise SystemExit('OU-II one-line enableTuner definition not found')
p.write_text(text.replace(old, new, 1))
