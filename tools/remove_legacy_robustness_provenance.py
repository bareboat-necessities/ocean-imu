#!/usr/bin/env python3
from pathlib import Path

path = Path(__file__).resolve().parent / "ou_robustness.py"
text = path.read_text(encoding="utf-8")
start = text.index("\ndef _restated_source_files(")
end = text.index("\n\ndef restat_bundle(", start)
text = text[:start] + text[end:]
path.write_text(text, encoding="utf-8")
print("removed obsolete moved-source restatement helper")
