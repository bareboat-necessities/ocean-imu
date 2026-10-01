import unittest
from tools.stability.ou3_theorem.outer_endpoint_bridge import audit
class T(unittest.TestCase):
 def test_no_invalid_separate_norm(self):
  r=audit(); self.assertEqual(r["causal_terminal_AW_reader_action_upper"],16); self.assertFalse(r["separate_endpoint_norm_product_valid_for_decisive_margin"]); self.assertIsNone(r["numeric_E_endpoint_plus_nonlinear"]); self.assertFalse(r["strict_margin_closed"])
if __name__=="__main__": unittest.main()
