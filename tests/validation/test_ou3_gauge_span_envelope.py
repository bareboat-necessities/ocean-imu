import math, unittest
from tools.stability.ou3_theorem.gauge_span_envelope import certificate,gauge_span_envelope
class GaugeSpanEnvelopeTests(unittest.TestCase):
 def test_current_constants(self):
  c=certificate()
  self.assertAlmostEqual(c["A_g_rad"],0.022961108159966083,places=14)
  self.assertAlmostEqual(c["peak_to_peak_cap_rad"],0.045922216319932166,places=14)
  self.assertAlmostEqual(c["L_g_rad_s"],0.001/9.80665,places=16)
  self.assertAlmostEqual(c["T_sat_s"],450.3431026738627,places=10)
  self.assertEqual(c["active_speed_bound"],"D_a/g")
 def test_branches(self):
  c=certificate(); self.assertAlmostEqual(gauge_span_envelope(100),c["L_g_rad_s"]*100)
  self.assertAlmostEqual(gauge_span_envelope(1000),c["peak_to_peak_cap_rad"])
if __name__=="__main__": unittest.main()
