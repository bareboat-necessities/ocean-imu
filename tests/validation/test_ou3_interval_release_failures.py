"""Expected enclosure failures must not hide implementation defects."""
import contextlib
import io
import json
import unittest
from unittest.mock import patch

from tools.stability.ou3_theorem import run_interval_release
from tools.stability.ou3_theorem.enclosure_failure import EnclosureFailure


class IntervalReleaseFailureTests(unittest.TestCase):
    def run_failure(self, error):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(run_interval_release, "root_cells", return_value=[object()]), \
                patch.object(run_interval_release, "propagate_history_cell", side_effect=error), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = run_interval_release.main()
        return status, json.loads(stdout.getvalue())

    def test_enclosure_failure_remains_non_promoting(self):
        status, result = self.run_failure(EnclosureFailure("guard", "interval straddles gate", .05))
        self.assertEqual(status, 2)
        self.assertEqual(result["classification"], "D_ENCLOSURE_FAILURE")
        self.assertEqual(result["first_failure_sample"], 10)
        self.assertFalse(result["verified"])

    def test_missing_history_field_is_an_implementation_failure(self):
        status, result = self.run_failure(KeyError("gravity_dir_0"))
        self.assertEqual(status, 1)
        self.assertEqual(result["classification"], "E_IMPLEMENTATION_FAILURE")
        self.assertIsNone(result["first_failure_sample"])
        self.assertFalse(result["verified"])


if __name__ == "__main__":
    unittest.main()
