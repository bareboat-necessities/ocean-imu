#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_brmm_physical_acceleration_witness_sector as S
class PhysicalAccelerationWitnessSectorTest(unittest.TestCase):
    def test_sector_closes_without_radial_conflation(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertTrue(d['joint_15D_source_quadratic_outer_relation_closed'])
        self.assertTrue(d['physical_source_scale_independent_of_hard_entry_radial'])
        self.assertTrue(d['zero_initial_error_nonzero_source_allowed'])
        self.assertFalse(d['hard_entry_radial_used_in_source_sector'])
        self.assertFalse(d['P4_PASS'])
    def test_canonical_sector_shapes(self):
        A0,A1,M,C=S.canonical_maps(.005);rows=S.physical_witness_sectors(A0,A1,M,C,8.0)
        self.assertEqual([n for n,_ in rows],['physical_a0_ball','physical_a1_ball','physical_J012_moment_iqc'])
        self.assertTrue(all((len(P),len(P[0]))==(16,16) for _,P in rows))
    def test_bad_detached_dimensions_rejected(self):
        z=S.I(0)
        with self.assertRaises(ValueError):S.physical_witness_sectors([[z]],[[z]],[[z]],[[z]],8.0)
if __name__=='__main__':unittest.main()
