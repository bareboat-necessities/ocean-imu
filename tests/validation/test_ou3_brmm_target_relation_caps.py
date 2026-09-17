"""The target relation must use the canonical cap, not historical literals."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/stability'))
import ou3_p4_complete_brmm_universal_target_relation as TARGET


class TargetRelationCapsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = TARGET.build()

    def test_both_caps_match_the_canonical_contract(self):
        cap = json.loads(TARGET.DEFAULT_DOMAIN.read_text())['normal_live'][
            'non_gravitational_cog_acceleration_norm_upper_mps2']
        self.assertEqual(self.report['declared_Normal_Live_acceleration_cap_mps2'], cap)
        self.assertEqual(self.report['qualified_source_acceleration_cap_mps2'], cap)
        self.assertEqual(TARGET.validate(self.report), [])

    def test_stale_or_nonfinite_caps_cannot_be_validated(self):
        for key in ('declared_Normal_Live_acceleration_cap_mps2',
                    'qualified_source_acceleration_cap_mps2'):
            for value in (4.0, 8.0, TARGET.ACCELERATION_CAP + 1, float('nan'), float('inf')):
                with self.subTest(key=key, value=value):
                    report = copy.deepcopy(self.report)
                    report[key] = value
                    self.assertTrue(TARGET.validate(report))

    def test_matching_but_noncanonical_caps_cannot_change_the_domain(self):
        report = copy.deepcopy(self.report)
        report['declared_Normal_Live_acceleration_cap_mps2'] = 4.0
        report['qualified_source_acceleration_cap_mps2'] = 4.0
        self.assertEqual(len(TARGET.validate(report)), 2)

    def test_cap_alignment_does_not_promote_P4(self):
        self.assertFalse(self.report['P4_promoted_here'])
        self.assertFalse(self.report['complete_source_cover_closed_here'])
        report = copy.deepcopy(self.report)
        report['P4_promoted_here'] = True
        self.assertIn('P4_promoted_here not false', TARGET.validate(report))


if __name__ == '__main__':
    unittest.main()
