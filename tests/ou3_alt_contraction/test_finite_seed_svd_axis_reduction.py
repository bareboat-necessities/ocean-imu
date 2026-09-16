from pathlib import Path
import os
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_seed_svd_axis_reduction as X


class Tests(unittest.TestCase):
    def test_every_post_QR_write_excludes_the_consumed_column(self):
        r = X.build()
        self.assertEqual(X.validate(r), [])
        self.assertEqual(r['all_Jacobi_rotation_column_pairs'], ((1,0),))
        self.assertEqual(r['all_sorting_swap_column_pairs'], ((0,1),))
        self.assertEqual(r['post_QR_V_write_columns'], (0,1))
        self.assertEqual(r['QR_householder_count'], 2)
        self.assertTrue(r['returning_solver_axis_equals_QR_Q_column_2'])
        self.assertTrue(r['rank_1_and_rank_2_same_column_preservation'])
        self.assertFalse(r['Jacobi_rotation_roundoff_accumulates_in_axis_column'])

    def test_returning_execution_does_not_prove_termination_or_axis_bound(self):
        r = X.build()
        for key in ('source_uniform_QR_axis_totality_and_bounds_closed',
                    'source_uniform_Jacobi_loop_termination_closed',
                    'target_compiler_arithmetic_correspondence_closed',
                    'near_antiparallel_JacobiSVD_solver_correspondence_qualified',
                    'storage_search_allowed'):
            self.assertFalse(r[key])
        r['source_uniform_Jacobi_loop_termination_closed'] = True
        self.assertTrue(X.validate(r))

    def test_source_header_identity_is_required(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError, 'unreviewed Eigen startup source'):
                X.audit_headers(td)
        # An explicit root allows a deployment source audit without making a
        # host Eigen install or a network fetch a prerequisite for unit tests.
        root = os.environ.get('OU3_ALT_TARGET_EIGEN_INCLUDE_DIR')
        if root:
            self.assertEqual(X.audit_headers(Path(root)), X.SOURCES)


if __name__ == '__main__':
    unittest.main()
