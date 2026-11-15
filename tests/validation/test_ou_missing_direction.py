"""Missing estimates survive serialization without becoming observations."""
import math
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import ou_evidence_provenance as provenance
import ou_validation as validation


class MissingDirectionTests(unittest.TestCase):
    def test_unresolved_sea_does_not_erase_other_direction_measurements(self):
        scenarios = ['stationary_jonswap_H0_270_case', 'stationary_jonswap_H1_500_case']
        summary = {
            (scenario, 'OU_III', 'Adaptive', 'dir_axis_abs_error_deg'): {'mean': value}
            for scenario, value in zip(scenarios, [math.nan, 12.0])
        }
        output = '\n'.join(validation._direction_table(scenarios, summary, lambda *a: '--'))
        self.assertIn(r'\OUValidationDirectionAbsError}{12.00}', output)
        self.assertIn(r'\OUValidationDirectionWorstAbsError}{12.00}', output)
        self.assertIn(r'\OUValidationTravelAbsError}{--}', output)
        self.assertNotIn('{nan}', output)

    def test_json_null_matches_csv_nan_but_never_a_measurement(self):
        self.assertTrue(provenance._scalar_matches_csv(None, 'nan'))
        self.assertTrue(provenance._scalar_matches_csv(None, ''))
        self.assertFalse(provenance._scalar_matches_csv(None, '0'))
        self.assertFalse(provenance._scalar_matches_csv(None, '12.5'))
        self.assertFalse(provenance._scalar_matches_csv(12.5, 'nan'))
        self.assertFalse(provenance._scalar_matches_csv(None, 'inf'))

    def test_summary_and_pairs_count_only_resolved_observations(self):
        left = [{'seed': 1, 'direction': 10.0}, {'seed': 2, 'direction': None},
                {'seed': 3, 'direction': math.nan}, {'seed': 4, 'direction': 20.0}]
        right = [{'seed': 1, 'direction': 8.0}, {'seed': 2, 'direction': 4.0},
                 {'seed': 3, 'direction': None}, {'seed': 4, 'direction': None}]
        np.testing.assert_array_equal(validation._finite_values(left, 'direction'), [10, 20])
        effect = validation._paired_effect(left, right, 'direction', ['seed'], 100,
                                           np.random.default_rng(7))
        self.assertEqual(effect['n_pairs'], 1)
        self.assertEqual(effect['mean_paired_difference'], 2.0)
        self.assertIsNone(validation._paired_effect(left[1:3], right, 'direction', ['seed'],
                                                   100, np.random.default_rng(7)))


if __name__ == '__main__':
    unittest.main()
