#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_source_indexed_phi_rebase as R

class SourceIndexedPhiRebaseTest(unittest.TestCase):
    def test_exact_rebase_algebra(self):
        d=R.build();self.assertEqual(R.validate(d),[])
        self.assertTrue(d['source_coordinate_rebase_exact_algebra_closed'])
        self.assertTrue(d['physical_state_unchanged_by_rebase'])
        self.assertTrue(d['C_equals_L_identity_for_rebase'])
        self.assertTrue(d['rebase_has_no_interior_epsilon_transport'])
        self.assertTrue(d['epsilon_packetwise_rezero_forbidden'])
        self.assertFalse(d['independent_rebase_disturbance_port_used'])
        self.assertFalse(d['source_uniform_rebase_defect_bound_closed_here'])
        self.assertFalse(d['production_rebases_inserted_between_all_samples_here'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])

if __name__=='__main__':unittest.main()
