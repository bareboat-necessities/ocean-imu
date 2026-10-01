import unittest
from tools.stability.ou3_theorem.a21_outer_entry import entry_reduction,conditional_entry
class T(unittest.TestCase):
 def test_fail_closed_until_physical_temporal_numbers_exist(self):
  r=entry_reduction(); self.assertTrue(r["compact_annulus_positive_homogeneous_loss_exists"]); self.assertFalse(r["finite_inner_entry_closed"])
 def test_conditional_margin(self):
  self.assertTrue(conditional_entry(loss_floor=2,supply_ceiling=1)["finite_entry"])
  self.assertFalse(conditional_entry(loss_floor=1,supply_ceiling=1)["finite_entry"])
if __name__=="__main__": unittest.main()
