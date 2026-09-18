import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.finite_error import (
    ServiceSuperwordWitness,
    audit_service_superword,
    incomplete_diagnostic,
)


class FiniteErrorTargetTests(unittest.TestCase):
    def test_point_diagnostic_remains_incomplete(self):
        w=incomplete_diagnostic("h",(1.0,0.95,0.90),1.0)
        r=audit_service_superword(w)
        self.assertTrue(r["strict_rho"])
        self.assertFalse(r["structural_premises"])
        self.assertFalse(r["magnetic_information"])
        self.assertFalse(r["certificate_complete"])
        self.assertFalse(r["point_witness_promoted_to_source_uniform_theorem"])

    def test_all_premises_are_required_for_complete_certificate(self):
        w=ServiceSuperwordWitness(
            history_id="h",
            V_start=1.0,
            V_end=0.8,
            disturbance_energy=0.0,
            rho_candidate=0.9,
            disturbance_gain=0.0,
            prefix_values=(1.0,0.9,0.8),
            prefix_bounds=(1.1,1.0,0.9),
            magnetic_information_min_eigenvalue=1.2,
            magnetic_information_floor=1.0,
            complete_shipping_map=True,
            inherited_state=True,
            same_history_motion=True,
            same_history_bias=True,
            applied_magnetic_information=True,
            finite_error_map=True,
            arithmetic_enclosed=True,
        )
        self.assertTrue(audit_service_superword(w)["certificate_complete"])

if __name__=="__main__":
    unittest.main()
