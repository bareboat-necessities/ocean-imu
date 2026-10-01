import math, unittest
from tools.stability.ou3_theorem.bias_gauge import certificate, qualified_fast_endpoint_envelope, slow_quiet_alias_envelope
class BiasGaugeTests(unittest.TestCase):
 def test_general_slow_only_bound(self):
  r=slow_quiet_alias_envelope(60,9.80665,.22516660498395405,.001,.02)
  self.assertAlmostEqual(r["joint_span_ceiling_rad"],2*math.asin(.06/(2*9.80665)))
  self.assertLess(r["joint_span_ceiling_rad"],math.radians(.351))
 def test_unknown_fast_stays_open(self):
  r=qualified_fast_endpoint_envelope(60,9.80665,.225,.001,.3)
  self.assertIsNone(r["span_ceiling_rad"]); self.assertEqual(r["temporal_qualification"],"OPEN")
 def test_qualified_fast_caps(self):
  r=qualified_fast_endpoint_envelope(60,9.80665,.225,.001,.3,.01,.02)
  self.assertAlmostEqual(r["span_ceiling_rad"],2*math.asin(.09/(2*9.80665)))
  with self.assertRaises(ValueError): qualified_fast_endpoint_envelope(60,9.80665,.225,.001,.3,.31,0)
 def test_fail_closed(self):
  r=certificate(); self.assertTrue(r["universal_slow_only_quiet_alias_bound_proved"])
  self.assertFalse(r["current_marine_slow_only_gauge_exclusion"]); self.assertFalse(r["full_slow_fast_joint_gauge_exclusion"])
if __name__=="__main__": unittest.main()
