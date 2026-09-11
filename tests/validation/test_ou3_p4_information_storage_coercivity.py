#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_information_storage_coercivity as C

class InformationStorageCoercivityTest(unittest.TestCase):
    def test_exact_boundary_compatibility_and_lower_coercivity_close(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        self.assertEqual(d['consecutive_word_metric_comparison_mu'],1.0)
        self.assertTrue(d['uniform_lower_coercivity_m_minus_closed'])
        for row in d['modes'].values():
            self.assertGreater(row['information_storage_coercivity_lower_m_minus'],0)
            self.assertGreater(row['prediction_Q_lambda_min_lower'],0)
    def test_upper_coercivity_stays_fail_closed_without_information_rate_bound(self):
        d=C.build();self.assertTrue(d['declared_async_magnetic_PE_has_no_uniform_accepted_information_rate_upper'])
        self.assertFalse(d['uniform_upper_coercivity_m_plus_closed']);self.assertFalse(d['uniform_storage_coercivity_closed'])
        self.assertFalse(d['hardware_magnetometer_ODR_substituted_as_theorem_premise'])
        self.assertIn('cumulative accepted-vector-information',d['weakest_missing_premise'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])
    def test_finite_information_budget_gives_strict_covariance_lower(self):
        for q in (1e-12,1e-8,1e-3):
            p=C.posterior_lower_from_information_budget(q,1000.0)
            self.assertGreater(p,0);self.assertLessEqual(p,q)

if __name__=='__main__':unittest.main()
