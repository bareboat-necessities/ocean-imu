from pathlib import Path

path = Path("doc/kalman_ou_iii/w3d-fus-methods.tex-part")
text = path.read_text(encoding="utf-8")
old = "% =======================================================\n\\section{Horizontal Wave Direction: Angle and Sign}"
new = "% =======================================================\n\\input{w3d-direction-methods.tex-part}\n"
if new in text:
    raise SystemExit(0)
if text.count(old) != 1:
    raise RuntimeError("expected obsolete wave-direction marker exactly once")
path.write_text(text.split(old, 1)[0] + new, encoding="utf-8")
