import unittest
from tools.stability.ou3_theorem.moving_gauge_obstruction import certificate
class MovingGaugeObstructionTests(unittest.TestCase):
    def test_exact_moving_pair_is_admitted(self):
        c=certificate()
        self.assertEqual(c["classification"],"B"); self.assertTrue(c["admitted_by_numeric_physical_envelopes"])
        self.assertTrue(c["delivered_accel_records_identical"] and c["delivered_gyro_records_identical"] and c["delivered_mag_records_identical"])
        self.assertTrue(c["fast_components_zero"]); self.assertLess(c["slow_accel_bias_norm_mps2"],0.22516660498395405)
        self.assertLess(c["slow_accel_bias_rate_max_mps3"],0.001)
        self.assertAlmostEqual(2*c["rock_amplitude_rad"],0.03490658503988659,places=15)
        self.assertAlmostEqual(2*c["position_amplitude_m"],0.03,places=15)
if __name__=="__main__": unittest.main()
