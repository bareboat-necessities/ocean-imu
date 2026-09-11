#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest,math
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_information_storage_coercivity as C

class InformationStorageCoercivityTest(unittest.TestCase):
    def test_exact_boundary_compatibility_and_uniform_coercivity_close(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        self.assertEqual(d['consecutive_word_metric_comparison_mu'],1.0)
        self.assertTrue(d['canonical_P3_following_prediction_boundary_consumed'])
        self.assertTrue(d['uniform_lower_coercivity_m_minus_closed'])
        self.assertTrue(d['uniform_upper_coercivity_m_plus_closed'])
        self.assertTrue(d['uniform_word_boundary_storage_coercivity_closed'])
        for row in d['modes'].values():
            self.assertGreater(row['information_storage_m_minus'],0)
            self.assertTrue(math.isfinite(row['information_storage_m_plus']))
            self.assertGreater(row['post_prediction_P_lambda_min_lower_from_Q'],0)

    def test_async_packet_rate_not_added_as_storage_premise(self):
        d=C.build()
        self.assertFalse(d['arbitrary_async_measurement_count_requires_boundary_m_plus_premise'])
        self.assertFalse(d['hardware_magnetometer_ODR_substituted_as_theorem_premise'])
        self.assertFalse(d['interior_Joseph_P_inverse_uniform_upper_bound_claimed'])
        self.assertFalse(d['word_boundary_moved_to_nonshipping_virtual_event'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])

    def test_process_lower_directly_bounds_information_metric_above(self):
        d=C.build()
        for row in d['modes'].values():
            q=row['post_prediction_P_lambda_min_lower_from_Q']
            self.assertGreater(q,0)
            self.assertGreaterEqual(row['information_storage_m_plus'],1.0/q)

if __name__=='__main__':unittest.main()
