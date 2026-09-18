from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[2]
LEGACY_TOKEN_RE=re.compile(r"\b(?:P[2-5]|BIAS[0-2]|BRMM)\b")
LEGACY_PATH_FRAGMENTS=["ou3_"+"p"+str(i) for i in range(2,6)]+["ou3_"+"br"+"mm"]+["bias"+str(i) for i in range(3)]
SKIP_TOP={".git","third_party","sim-data-files"}
BINARY_SUFFIXES={".pdf",".png",".jpg",".jpeg",".gif",".svg",".xz",".gz",".zip",".bin",".ico",".woff",".woff2",".ttf"}

class ArchitectureCleanupTests(unittest.TestCase):
    def test_no_retired_architecture_survives_repository(self):
        bad=[]
        for path in ROOT.rglob("*"):
            if not path.is_file():
                continue
            rel_path=path.relative_to(ROOT)
            if rel_path.parts and rel_path.parts[0] in SKIP_TOP:
                continue
            rel=rel_path.as_posix()
            rel_lower=rel.lower()
            if any(x in rel_lower for x in LEGACY_PATH_FRAGMENTS):
                bad.append("path:"+rel)
                continue
            if path.suffix.lower() in BINARY_SUFFIXES:
                continue
            try:
                text=path.read_text(encoding="utf-8")
            except (UnicodeDecodeError,OSError):
                continue
            match=LEGACY_TOKEN_RE.search(text)
            if match:
                bad.append(f"text:{rel}:{match.group(0)}")
        self.assertEqual([],bad,"retired proof architecture survived:\n"+"\n".join(bad[:200]))

    def test_exact_old_proof_workflow_names_are_absent(self):
        workflows=ROOT/".github/workflows"
        old=[
            p.name for p in workflows.glob("ou3-*.yml")
            if p.name not in {"ou3-stability-proof.yml","ou3-lever-arm-study.yml"}
        ]
        self.assertEqual([],sorted(old))

if __name__=="__main__":
    unittest.main()
