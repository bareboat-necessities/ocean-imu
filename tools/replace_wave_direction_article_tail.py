from pathlib import Path

path = Path("doc/kalman_ou_iii/w3d-fus-methods.tex-part")
text = path.read_text(encoding="utf-8")
marker = "% =======================================================\n\\section{Horizontal Wave Direction: Angle and Sign}"
if text.count(marker) != 1:
    raise RuntimeError("expected obsolete wave-direction marker exactly once")
prefix = text.split(marker, 1)[0]
path.write_text(
    prefix + "% =======================================================\n"
    + "\\input{w3d-direction-methods.tex-part}\n",
    encoding="utf-8",
)
