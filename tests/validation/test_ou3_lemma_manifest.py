import json,re,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
MAN=ROOT/"reports/results/ou3_stability/lemma-manifest.json"
TEX=ROOT/"doc/kalman_ou_iii/kalman_ou-w3d-stability-study.tex"
class LemmaManifestTests(unittest.TestCase):
 def load(self):return json.loads(MAN.read_text())
 def test_exact_latex_coverage(self):
  m=self.load();text=TEX.read_text()
  found=[(k,t) for k,t in re.findall(r"\\begin\{(lemma|theorem)\}\[([^\]]+)\]",text)]
  declared=[(x["kind"],x["latex_title"]) for x in m["lemmas"]]
  self.assertEqual(found,declared)
  self.assertEqual(len({x["id"] for x in m["lemmas"]}),len(declared))
 def test_status_contract(self):
  for x in self.load()["lemmas"]:
   self.assertIn(x["status"],("PROVED","CONDITIONAL","OPEN"))
   if x["status"]=="CONDITIONAL":self.assertTrue(x["open_dependencies"])
   if x["status"]=="OPEN":self.assertTrue(x["open_dependencies"])
   self.assertTrue(x["assumptions"]);self.assertTrue(x["proof_modules"])
 def test_referenced_paths_exist(self):
  for x in self.load()["lemmas"]:
   for p in x["proof_modules"]+x["certificates"]:
    self.assertTrue((ROOT/p).exists(),f"{x['id']}: missing {p}")
if __name__=="__main__":unittest.main()
