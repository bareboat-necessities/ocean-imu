"""The actual entry point must not start a coarse search on a Jacobian word."""
import unittest
from unittest.mock import patch
from tools.stability.ou3_alt_contraction import first_common_metric_attempt as M


class FirstCommonMetricAttemptTest(unittest.TestCase):
    def test_candidate_helper_cannot_bypass_the_finite_master_guard(self):
        with self.assertRaisesRegex(RuntimeError, 'finite-state storage blocked'):
            M.normalized_diagonal_metric()

    def test_first_common_metric_attempt_requires_a_finite_master(self):
        with patch.object(M, 'normalized_diagonal_metric') as candidate:
            with self.assertRaisesRegex(RuntimeError, 'finite-state storage blocked'):
                M.build()
            candidate.assert_not_called()


if __name__ == '__main__':
    unittest.main()
