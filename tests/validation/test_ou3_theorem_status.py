import json,sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.theorem_status import status_report

class TheoremStatusTests(unittest.TestCase):
    def test_single_architecture_is_fail_closed(self):
        r=status_report()
        self.assertEqual(r["principal_assumptions"],["MARINE MOTION","IMU BIAS","MAGNETIC SERVICE"])
        self.assertFalse(r["parallel_no_magnetometer_stability_path"])
        self.assertFalse(r["theorem_closed"])
        self.assertFalse(r["regional_practical_stability_claimed"])
        self.assertIn("H18-to-A21 release",r["proof_path"])
    def test_committed_status_matches_code(self):
        p=ROOT/"reports/results/ou3_stability/theorem-status.json"
        self.assertEqual(json.loads(p.read_text()),status_report())

if __name__=="__main__": unittest.main()
