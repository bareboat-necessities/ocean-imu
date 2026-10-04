import struct,unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.planar_service_frontend_binding import (
    dyadic,round32,inverse_sqrt_word,certificate,ShippingBindingFailure)
from tools.stability.ou3_theorem.planar_service_guard import certificate as guard
from tools.stability.ou3_theorem.planar_service_mahony_tube import certificate as tube
from tools.stability.ou3_theorem.causal_tuner_interval import I
from tools.stability.ou3_theorem.literal_history_leaf_propagator import propagate_history_cell

class PlanarServiceFrontendTests(unittest.TestCase):
    def test_binary32_exact_rounding_and_ties(self):
        for x in [F(1),F(1,200),F(-1,3),F(1,2**149),F(1,2**150),F(196133,20000)]:
            expected=struct.unpack('<I',struct.pack('<f',float(x)))[0]
            self.assertEqual(round32(x),expected)
        self.assertEqual(round32(F(1)+F(1,2**24)),0x3f800000)
        self.assertEqual(round32(F(1)+F(3,2**24)),0x3f800002)

    def test_fast_normalization_is_not_exact_unit_normalization(self):
        bits,steps=inverse_sqrt_word(0x3f800000)
        self.assertLess(dyadic(bits)**2,F(1))
        self.assertGreater(dyadic(bits),F(99,100))
        self.assertEqual(len(steps),7)
        c=certificate()
        self.assertFalse(c['unit_norm_after_shipping_float_normalization'])
        self.assertTrue(c['pre_normalization_quaternion_identity_retained'])
        self.assertFalse(c['shipping_frontend_interval_binding_verified'])

    def test_reference_frontend_cannot_generate_shipping_leaf(self):
        with self.assertRaises(ShippingBindingFailure):propagate_history_cell(None)

    def test_interval_clamp_handles_entirely_outside_boxes(self):
        self.assertEqual(I(-3,-1).clamp(0,1),I(0,0))
        self.assertEqual(I(2,4).clamp(0,1),I(1,1))
        self.assertEqual(I(-2,4).clamp(0,1),I(0,1))
        with self.assertRaises(ValueError):I(0,1).clamp(2,1)

    def test_planar_guard_all_time_real_arithmetic_margin(self):
        c=guard()
        self.assertLess(F(c['second_detector_output_norm_upper']),F(c['engage_floor_mps2']))
        self.assertGreater(F(c['exp_lower_polynomial']),1/F(c['detector_decay_upper']))
        self.assertLess(c['detector_RMS_upper_mps2'],.000531)
        self.assertFalse(c['float32_guard_transfer_verified'])
        self.assertFalse(c['all_time_magnetic_service_verified'])

    def test_conditional_mahony_tube_keeps_error_binding_open(self):
        c=tube()
        self.assertGreater(F(c['strict_radial_margin']),0)
        for row in c['endpoint_checks']:
            self.assertGreater(min(map(F,row['positive_LDL_pivots'])),0)
        self.assertFalse(c['complete_literal_one_step_error_caps_verified'])
        self.assertFalse(c['all_time_frontend_tube_verified'])
        self.assertTrue(c['raw_quaternion_norm_not_assumed_one'])

if __name__=='__main__':unittest.main()
