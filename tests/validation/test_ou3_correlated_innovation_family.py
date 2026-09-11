from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"tools/stability"))

from ou3_interval import Interval, matrix_point
import ou3_correlated_innovation_family as C


class CorrelatedInnovationFamilyTests(unittest.TestCase):
    def test_same_triplet_owns_S_inverse_and_K(self):
        P=matrix_point([[2.0,.2,0],[.2,1.5,.1],[0,.1,1.0]])
        H=matrix_point([[1,0,0],[0,1,0],[0,0,1]])
        R=matrix_point([[.2,0,0],[0,.3,0],[0,0,.4]])
        f=C.build(P,H,R)
        p=f.inverse_provenance
        self.assertEqual(p["primitive_family"], "(P,H,R)")
        self.assertTrue(p["same_P_H_R_used_for_PHt_S_Sinv_K"])
        self.assertFalse(p["S_is_external_independent_argument"])
        self.assertFalse(p["rectangular_S_family_can_be_selected_independently"])
        self.assertEqual(p["K_identity"], "K=P H^T (H P H^T+R)^-1")
        self.assertFalse(p["P4_PASS"])

    def test_wide_S_hull_does_not_become_an_external_family(self):
        # H uncertainty can make an entrywise S hull very wide.  The API still
        # takes only P,H,R and derives S internally; callers cannot substitute a
        # singular rectangular S member unrelated to the same triplet.
        P=matrix_point([[1,0,0],[0,1,0],[0,0,1]])
        z=Interval(-1.0,1.0)
        H=[[z,z,z],[z,z,z],[z,z,z]]
        R=matrix_point([[.1,0,0],[0,.1,0],[0,0,.1]])
        f=C.build(P,H,R)
        self.assertTrue(f.inverse_provenance["same_P_H_R_used_for_PHt_S_Sinv_K"])
        self.assertFalse(f.inverse_provenance["S_is_external_independent_argument"])


if __name__ == "__main__":
    unittest.main()
