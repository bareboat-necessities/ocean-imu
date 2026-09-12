from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_private_mahony_discrete_comparison as mod  # noqa: E402
import ou3_brmm_private_mahony_live_invariant as CONT  # noqa: E402


class BrmmPrivateMahonyDiscreteComparisonTest(unittest.TestCase):
    """The discrete comparison is blocked by its continuous prerequisite.

    The padded 8.8 m/s^2 envelope leaves the 87 deg-chart PI invariant with an
    empty admissible level window, so the discrete/binary32 comparison refuses
    to build rather than composing a step on top of an unclosed invariant.  The
    refusal is the certificate here; the discrete lemma cannot be reinstated
    until the continuous invariant closes.
    """

    def test_refuses_to_compose_on_an_unclosed_continuous_invariant(self):
        failures = CONT.validate(CONT.build())
        self.assertNotEqual(failures, [])
        with self.assertRaisesRegex(RuntimeError, "continuous Mahony invariant invalid"):
            mod.build()

    def test_refusal_names_the_seed_level_window_obstruction(self):
        with self.assertRaises(RuntimeError) as ctx:
            mod.build()
        message = str(ctx.exception)
        self.assertIn("continuous_all_live_PI_invariant_closed is not true", message)
        self.assertIn("no metric level contains the seed", message)


if __name__ == "__main__":
    unittest.main()
