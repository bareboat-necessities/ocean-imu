"""Speed comparisons must not silently change executables or omit a sea."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import encounter_regression_inventory as inventory


class EncounterInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = ROOT / 'reports/results/filter_sim_gates'
        cls.moving = json.loads((folder / 'measurements.json').read_text())
        cls.stationary = json.loads((folder / 'stationary_control.json').read_text())

    def test_pairs_by_record_identity_not_report_order(self):
        expected = inventory.compare(self.moving, self.stationary)
        control = copy.deepcopy(self.stationary)
        for family in inventory.FAMILIES:
            control['families'][family]['records'].reverse()
        self.assertEqual(inventory.compare(self.moving, control), expected)
        self.assertEqual(len(expected['rows']), 64)
        row = expected['rows'][0]['metrics']['disp_z_pct_hs']
        self.assertAlmostEqual(row['delta'], row['moving'] - row['stationary'])

    def test_rejects_changed_binary_incomplete_seas_and_invalid_metrics(self):
        family = inventory.FAMILIES[0]
        for mutation in ('binary', 'missing', 'duplicate', 'nan', 'archive'):
            control = copy.deepcopy(self.stationary)
            block = control['families'][family]
            if mutation == 'binary':
                block['binary_sha256'] = 'different executable'
            elif mutation == 'missing':
                block['records'].pop()
            elif mutation == 'duplicate':
                block['records'][0] = block['records'][1]
            elif mutation == 'nan':
                block['records'][0]['disp_z_pct_hs'] = float('nan')
            else:
                control['archive_sha256'] = 'obsolete dataset'
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                inventory.compare(self.moving, control)


if __name__ == '__main__':
    unittest.main()
