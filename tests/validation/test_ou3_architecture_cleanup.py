from __future__ import annotations

from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[2]
UNIQUE_LEGACY=["BR"+"MM"]+["B"+"IAS"+str(i) for i in range(3)]
STAGE_LEGACY=["P"+str(i) for i in range(2,6)]
LEGACY_PATH_FRAGMENTS=(
    ["ou3_"+"p"+str(i) for i in range(2,6)]
    +["ou3-"+"p"+str(i) for i in range(2,6)]
    +["ou3_"+"br"+"mm","ou3-"+"br"+"mm"]
    +["bias"+str(i) for i in range(3)]
)
SKIP_TOP={".git","third_party","sim-data-files"}
BINARY_SUFFIXES={".pdf",".png",".jpg",".jpeg",".gif",".svg",".xz",".gz",".zip",".bin",".ico",".woff",".woff2",".ttf"}
TEXT_CONTEXT_MARKERS=("ou3","kalman_ou_iii","kalman-ou-w3d-stability","stability")

def standalone(term: str,text: str) -> bool:
    return re.search(
        r"(?<![A-Za-z0-9])"+re.escape(term)+r"(?![A-Za-z0-9])",
        text,
        flags=re.IGNORECASE,
    ) is not None

class ArchitectureCleanupTests(unittest.TestCase):
    def test_no_retired_architecture_survives_repository(self):
        bad=[]
        for path in ROOT.rglob("*"):
            if not path.is_file(): continue
            rel_path=path.relative_to(ROOT)
            if rel_path.parts and rel_path.parts[0] in SKIP_TOP: continue
            rel=rel_path.as_posix(); low=rel.lower()
            if any(x in low for x in LEGACY_PATH_FRAGMENTS):
                bad.append("path:"+rel); continue
            if path.suffix.lower() in BINARY_SUFFIXES: continue
            try: text=path.read_text(encoding="utf-8")
            except (UnicodeDecodeError,OSError): continue
            for term in UNIQUE_LEGACY:
                if standalone(term,text): bad.append(f"text:{rel}:{term}")
            proof_context=(
                rel=="AGENTS.md"
                or low.startswith(".github/workflows/")
                or any(marker in low for marker in TEXT_CONTEXT_MARKERS)
                or "OU-III" in text
                or "OU3" in text
            )
            if proof_context:
                for term in STAGE_LEGACY:
                    if standalone(term,text): bad.append(f"text:{rel}:{term}")
        self.assertEqual([],sorted(set(bad)),
                         "retired OU-III stability architecture survived:\n"
                         +"\n".join(sorted(set(bad))[:200]))

    def test_exact_old_proof_workflow_names_are_absent(self):
        workflows=ROOT/".github/workflows"
        old=[
            p.name for p in workflows.glob("ou3-*.yml")
            if p.name not in {"ou3-stability-proof.yml","ou3-lever-arm-study.yml"}
        ]
        self.assertEqual([],sorted(old))

if __name__=="__main__":
    unittest.main()
