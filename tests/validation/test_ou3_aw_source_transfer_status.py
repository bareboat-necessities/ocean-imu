import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.aw_source_transfer_status import certificate

class T(unittest.TestCase):
    def test_conservative_real_source_budget_has_large_Astar_room(self):
        r=certificate()
        self.assertGreater(F(r["remaining_Astar_margin_before_float32_mps2"]),F("0.881"))
        self.assertFalse(r["real_arithmetic_reset_error_dependent_twist_in_source"])
        self.assertTrue(r["Astar_real_arithmetic_source_budget_closed"])
        self.assertTrue(r["Astar_source_uniform_strictness_closed"])
        self.assertFalse(r["quaternion_polynomial_is_separate_AW_source"])
        self.assertFalse(r["explicit_G0_point4_point12_premises_proved"])
        # The source theorem is deliberately weaker than the old explicit G0
        # constants but is strict enough to exclude persistent field-compatible
        # zero-dissipation geometry on qualified moving superwords.

if __name__=="__main__":
    unittest.main()
