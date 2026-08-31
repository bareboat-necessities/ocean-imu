from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_outer_sector_capture_certificate as P5

WORKFLOW = ROOT / ".github" / "workflows" / "ou3-usable-certificates-fast.yml"


class Ou3CertificateDeadEndRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sector = SECTOR.build()
        cls.p5 = P5.build()
        cls.workflow = WORKFLOW.read_text(encoding="utf-8")

    def test_microscopic_sector_can_never_be_promoted_as_usable(self):
        d = deepcopy(self.sector)
        d["design_full_attitude_angle_rad"] = 0.00124
        self.assertNotEqual(SECTOR.validate(d), [])

    def test_retired_global_scalar_accounting_flags_fail_closed(self):
        for key in (
            "global_packet_count_times_lipschitz_defect_used",
            "whole_word_weakest_P3_delta_used_as_attitude_sector_margin",
        ):
            with self.subTest(key=key):
                d = deepcopy(self.sector)
                d[key] = True
                self.assertNotEqual(SECTOR.validate(d), [])

        for key in (
            "legacy_uniform_transport_route_used",
            "legacy_microscopic_inner_seed_used_as_outer_capture_target",
        ):
            with self.subTest(key=key):
                d = deepcopy(self.p5)
                d[key] = True
                self.assertNotEqual(P5.validate(d), [])

    def test_partial_linear_improvement_cannot_claim_complete_p4(self):
        d = deepcopy(self.sector)
        d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"] = True
        self.assertNotEqual(SECTOR.validate(d), [])
        d = deepcopy(self.p5)
        d["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"] = True
        self.assertNotEqual(P5.validate(d), [])

    def test_pr441_route_ceiling_regression_remains_in_focused_ci(self):
        self.assertIn("test_ou3_p4_p5_route_ceiling", self.workflow)
        self.assertIn("tools/ou3_p4_p5_route_ceiling_certificate.py", self.workflow)

    def test_useful_438_full_word_progress_is_wired_but_not_promoted(self):
        self.assertIn("tools/ou3_p4_translation_full_word_rigorous.py", self.workflow)
        self.assertIn("tools/ou3_p4_post_translation_bottleneck.py", self.workflow)
        self.assertIn("test_ou3_p4_post_translation_bottleneck", self.workflow)

    def test_retired_438_scalar_frontiers_are_not_promotion_inputs(self):
        # These files may remain in old branches/history as diagnostics, but the
        # usable-certificate workflow must never gate theorem promotion on them.
        forbidden = (
            "ou3_p4_direct_word_contraction_certificate.py",
            "ou3_p4_nextgen_gain_certificate.py",
            "ou3_p4_nextgen_directional_certificate.py",
            "ou3_p4_nextgen_widened_certificate.py",
            "ou3_p4_thirdgen_combined_certificate.py",
            "ou3_p4_frontier_combined_certificate.py",
        )
        for name in forbidden:
            with self.subTest(name=name):
                self.assertNotIn(name, self.workflow)


if __name__ == "__main__":
    unittest.main()
