import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2] / 'tools' / 'stability'


class Tests(unittest.TestCase):
    def _read(self, name):
        return json.loads((ROOT / name).read_text())

    def test_fcr_audit_remains_an_explicit_execution_premise(self):
        report = self._read('ou3_alt_target_fcr_runtime.json')
        self.assertEqual(report['architectural_FCR_reset_value'], 'undefined')
        self.assertTrue(report['task_first_FPU_use_inherits_active_FCR'])
        self.assertFalse(report['rounding_mode_runtime_initialization_and_preservation_proved'])
        self.assertIn('both cores enter the admitted execution with FCR.RM=0',
                      report['required_execution_premises'])

    def test_math_link_and_pow_range_are_not_overpromoted(self):
        link = self._read('ou3_alt_target_math_link.json')
        powf = self._read('ou3_alt_target_powf_range.json')
        self.assertTrue(link['standalone_exp_log_sqrt_namespace_resolution_qualified'])
        self.assertFalse(link['whole_firmware_link_resolution_qualified'])
        self.assertTrue(powf['range_certificate']['all_intermediates_finite'])
        self.assertTrue(powf['range_certificate']['return_strictly_positive'])
        self.assertFalse(powf['range_certificate']['approximation_accuracy_qualified'])


if __name__ == '__main__':
    unittest.main()
