import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.execution_contract import audit_same_execution
from tools.stability.ou3_theorem.imu_bias import BiasSample
from tools.stability.ou3_theorem.magnetic_service import MagneticEvent
from tools.stability.ou3_theorem.marine_motion import MarineSample


class SameExecutionTests(unittest.TestCase):
    def marine(self, history):
        return [
            MarineSample(
                history,
                0.0,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
            )
        ]

    def bias(self, history):
        return [BiasSample(history, 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))]

    def mag(self, history):
        return [
            MagneticEvent(
                history,
                0.0,
                True,
                True,
                True,
                False,
                ((1.0, 0.0), (0.0, 1.0)),
            )
        ]

    def test_all_principal_contracts_share_one_execution(self):
        r = audit_same_execution(self.marine("h"), self.bias("h"), self.mag("h"))
        self.assertTrue(r["same_execution_pass"], r["failures"])
        self.assertEqual(r["execution_id"], "h")
        self.assertFalse(r["independent_per_word_reselection_allowed"])

    def test_detached_magnetic_history_fails(self):
        r = audit_same_execution(self.marine("h"), self.bias("h"), self.mag("other"))
        self.assertFalse(r["same_execution_pass"])

    def test_missing_principal_contract_fails(self):
        r = audit_same_execution(self.marine("h"), self.bias("h"), [])
        self.assertFalse(r["same_execution_pass"])
        self.assertIn("magnetic_service: missing principal-assumption records", r["failures"])


if __name__ == "__main__":
    unittest.main()
