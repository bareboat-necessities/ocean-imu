from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
SCOPES=[ROOT/"AGENTS.md",ROOT/".github/workflows",ROOT/"doc/kalman_ou_iii",ROOT/"docs",
        ROOT/"tools/stability",ROOT/"tests/validation",ROOT/"reports/results/ou3_stability"]
LEGACY_TERMS=["P"+str(i) for i in range(2,6)]+["B"+"IAS"+str(i) for i in range(3)]+["BR"+"MM"]
LEGACY_PATH_FRAGMENTS=["ou3_"+"p"+str(i) for i in range(2,6)]+["ou3_"+"br"+"mm"]+["bias"+str(i) for i in range(3)]

class ArchitectureCleanupTests(unittest.TestCase):
    def test_no_retired_architecture_survives_proof_surfaces(self):
        bad=[]
        for scope in SCOPES:
            paths=[scope] if scope.is_file() else list(scope.rglob("*"))
            for path in paths:
                if not path.is_file() or path.name==Path(__file__).name: continue
                rel=path.relative_to(ROOT).as_posix().lower()
                if any(x in rel for x in LEGACY_PATH_FRAGMENTS):
                    bad.append("path:"+rel); continue
                try: text=path.read_text(encoding="utf-8")
                except UnicodeDecodeError: continue
                for term in LEGACY_TERMS:
                    if term in text: bad.append(f"text:{rel}:{term}")
        self.assertEqual([],bad,"retired proof architecture survived:\n"+"\n".join(bad[:100]))

if __name__=="__main__": unittest.main()
