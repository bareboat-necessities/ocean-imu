#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_source_indexed_prefix_transport as T

class H18SourceIndexedPrefixTransportTest(unittest.TestCase):
    def test_exact_two_sample_transport(self):
        d=T.build();self.assertEqual(T.validate(d),[])
        self.assertTrue(d['H18_every_prefix_transport_constructor_closed'])
        self.assertTrue(d['cross_sample_Phi_rebase_consumed'])
        self.assertTrue(d['prediction_uses_trusted_shipping_F_not_nonlinear_Jacobian'])
        self.assertTrue(d['Joseph_C_equals_G_I_minus_KH'])
        self.assertTrue(d['Joseph_L_equals_G'])
        self.assertTrue(d['only_accelerometer_has_interior_epsilon_transport'])
        self.assertFalse(d['H18_endpoint_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])

if __name__=='__main__':unittest.main()
