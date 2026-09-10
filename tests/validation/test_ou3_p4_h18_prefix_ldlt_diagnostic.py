#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_prefix_ldlt_diagnostic as D
class H18PrefixLdltDiagnosticTest(unittest.TestCase):
    def test_diagnostic_uses_actual_prefix_matrices_without_promotion(self):
        d=D.build();self.assertEqual(D.validate(d),[])
        self.assertTrue(d['actual_synchronized_H18_prefix_transport_consumed'])
        self.assertTrue(d['actual_event_covariance_information_matrices_consumed'])
        self.assertTrue(d['full_augmented_interval_LDLT_backend_consumed'])
        self.assertTrue(d['one_conservative_independent_source_column_used_for_diagnostic'])
        self.assertFalse(d['source_moment_radial_correlation_preserved_in_diagnostic'])
        self.assertFalse(d['binary32_channel_included_in_diagnostic'])
        self.assertFalse(d['production_endpoint_LDLT_promoted_here'])
        self.assertFalse(d['P4_PASS'])
        self.assertGreater(len(d['prefixes']),0)
if __name__=='__main__':unittest.main()
