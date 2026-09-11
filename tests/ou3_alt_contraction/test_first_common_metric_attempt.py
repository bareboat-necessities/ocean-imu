import unittest
from tools.stability.ou3_alt_contraction import first_common_metric_attempt as M


class FirstCommonMetricAttemptTest(unittest.TestCase):
    def test_first_common_metric_attempt_is_fail_closed_and_universal(self):
        d=M.build();self.assertEqual(M.validate(d),[])
        print('ALT_FIRST_COMMON_METRIC',M.summary(d))
        self.assertTrue(d['same_M_used_for_all_H_A_BIAS_families'])
        self.assertFalse(d['trajectory_or_replay_fit_used_for_M'])
        self.assertFalse(d['per_mode_metric_used'])
        self.assertFalse(d['per_bias_family_metric_used'])
        self.assertFalse(d['async_magnetometer_nonexpansive_same_M_closed'])
        self.assertFalse(d['common_joint24_storage_complete'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':unittest.main()
