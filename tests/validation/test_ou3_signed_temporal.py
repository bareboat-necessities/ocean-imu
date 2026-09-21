import sys,unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.signed_temporal import *

class SignedTemporalTests(unittest.TestCase):
    def test_spline_atoms_and_endpoint_jets(self):
        k=(F(0),F(1),F(2),F(3))
        self.assertEqual(sum(atoms(k)),0)
        self.assertTrue(endpoint_jets_vanish(k))
    def test_homogeneous_adjoint_failure_is_exact(self):
        self.assertTrue(zero_terminal_homogeneous_adjoint_impossible(3))
    def test_physical_margin_positive(self):
        b=physical_bounds()
        self.assertGreater(b["physical_floor_margin"],0)
    def test_projection_gap_positive(self):
        self.assertGreater(projection_sector_gap(),0)
    def test_fail_closed(self):
        c=certificate()
        self.assertFalse(c["source_uniform_nominal_force_field_temporal_margin"])
        self.assertFalse(c["source_uniform_nominal_gyro_alias_temporal_margin"])
        self.assertFalse(c["uniform_historical_AG_readout_action"])
        self.assertFalse(c["full_21_covariance_upper"])
        self.assertFalse(c["rho0_certified"])

if __name__=="__main__": unittest.main()
