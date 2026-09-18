import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.no_mag_obstruction import obstruction_report,repeated_state

class NoHeadingServiceObstructionTests(unittest.TestCase):
    def test_unipotent_block_has_unit_spectral_radius(self):
        r=obstruction_report(3.0)
        self.assertEqual(r["spectral_radius"],1.0)
        self.assertFalse(r["strict_full_state_contraction_possible_without_heading_information"])
    def test_axial_bias_error_causes_linear_heading_growth(self):
        theta,bias=repeated_state(0.0,0.01,2.0,50)
        self.assertAlmostEqual(theta,1.0); self.assertAlmostEqual(bias,0.01)

if __name__=="__main__": unittest.main()
