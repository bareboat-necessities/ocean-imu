"""Report-logic tests for the consolidated P1--P5 certificate numbers.

The producer itself rebuilds every stage and takes minutes; the dedicated
``ou3-p1-p5-certificate-numbers`` workflow runs it.  What is checked here is the
part that decides what the report *says*: the usability thresholds, the verdict
rules, and the refusal to claim a complete proof while an obligation is open.
"""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p1_p5_certificate_numbers as REPORT


def _ctx():
    """A synthetic context carrying the numbers the shipping stages report."""
    return {
        "startup": {
            "operating_domain": {
                "latent_acceleration_error_fraction_g": 0.3,
                "physical_handoff_coordinate_bounds": {
                    "position_error_norm_upper_m": 20.0,
                    "velocity_error_norm_upper_mps": 5.0,
                    "latent_acceleration_error_norm_upper_mps2": 2.941995,
                    "gyro_bias_error_norm_upper_rad_s": 0.01,
                    "accelerometer_bias_error_norm_upper_mps2": 0.5,
                    "integral_displacement_error_norm_upper_m_s": 300.0,
                },
            },
        },
        "startup_validation_pass": True,
        "P1_normal_gauged_cayley_norm_upper": 0.2721648148683776,
        "P1_timeout_gauged_cayley_norm_upper": 0.5947333355555983,
        "P2": {
            "P2_SOURCE_PATH_CERTIFICATE": "PASS",
            "source_only": True,
            "trajectory_replay_used": False,
            "partition": {"states": 800},
            "recurrent_states": 800,
            "transition_edges": 640000,
            "strongly_connected_components": 1,
            "raw_tuner_sigma_subfloor_states_included": True,
            "RS_target_full_deployed_clamp_overapprox": True,
            "filter_sigma_floor_mps2": 0.05,
            "raw_tuner_sigma_partition_lower": 1e-06,
        },
        "word": {
            "modes": {
                "H": {"P3_word_endpoint_delta_lower": 2.2953997386276688e-20,
                      "P3_homogeneous_prefix_information_gain_upper": 1.0},
                "A": {"P3_word_endpoint_delta_lower": 2.2953997386276688e-20,
                      "P3_homogeneous_prefix_information_gain_upper": 1.0},
            },
        },
        "P4_sector": {
            "P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE": "PASS",
            "design_full_attitude_angle_rad": 0.8,
            "design_full_attitude_angle_deg": 45.836623610465864,
            "design_cayley_norm_upper": 0.845586437476324,
            "exact_vector_strong_monotonicity_factor_lower": 0.8483533546735822,
            "exact_eta_to_rotational_residual_information_ratio_upper": 0.17875410581097537,
            "P1_overlap": {
                "normal_gauged_cayley_norm_upper": 0.2721648148683776,
                "timeout_gauged_cayley_norm_upper": 0.5947333355555983,
                "normal_gauged_inside_sector": True,
                "timeout_gauged_inside_sector": True,
            },
        },
        "P4_first_accel_budget": {
            "declared_aw_over_lowest_specific_force": 0.5883990000000002,
            "ladder_rows": [{
                "angle_deg": 30.0, "ladder_family": "descending",
                "sector_invariance_correction_budget_upper_rad": 0.27225152012902093,
                "nuisance_correction_norm_upper_rad": 1.375114933176711,
                "nuisance_over_budget_ratio": 5.050899008846818,
            }],
            "smallest_probe_angle_deg": 1.0,
            "smallest_probe_nuisance_over_budget_ratio": 1.337951721023343,
            "shrinking_the_candidate_angle_alone_can_close_the_budget": False,
        },
        "P5_outer": {
            "P5_OUTER_SECTOR_CAPTURE_CERTIFICATE": "PASS",
            "declared_P5_entrance": {
                "gauged_full_attitude_angle_upper_deg": 45.0,
                "attitude_geometry": {"cayley_norm_upper": 0.8284271247461905},
                "position_component_abs_error_upper_Hs_factor": 0.5,
                "position_norm_upper_Hs_factor": 0.8660254037844388,
            },
            "outer_sector_angle_deg": 45.836623610465864,
            "all_source_handoff_branches_enter_outer_sector": True,
            "N_outer_words": 0,
            "P1_conservative_handoff_box_replaced": False,
            "upper_cosine_enclosure_used_for_ungauged_inclusion": True,
            "legacy_microscopic_inner_seed_used_as_outer_capture_target": False,
            "branches": {
                "normal_gauged": {"P1_handoff_cayley_norm_upper": 0.2721648148683776},
                "timeout_gauged": {"P1_handoff_cayley_norm_upper": 0.5947333355555983},
                "timeout_ungauged": {"full_heading_radius_assigned": False},
            },
        },
        "route_ceiling": {
            "modes": {m: {
                "certified_attitude_capture_radius_now": 5.808479596010723e-32,
                "route_ceiling_at_shipping_prefix_factor": 0.0012376237623762383,
                "route_can_reach_P1_handoff": False,
            } for m in ("H", "A")},
        },
        "translation": None,
        "first_accel": None,
    }


class StageVerdictTests(unittest.TestCase):
    def test_shipping_numbers_make_every_stage_geometry_usable(self):
        ctx = _ctx()
        stages = [REPORT._p1_stage(ctx), REPORT._p2_stage(ctx), REPORT._p3_stage(ctx),
                  REPORT._p4_stage(ctx), REPORT._p5_stage(ctx)]
        for s in stages:
            for c in s["usability_checks"]:
                self.assertTrue(c["pass"], f"{s['stage']}: {c['check']}")
        self.assertEqual([s["verdict"] for s in stages], [
            "USABLE", "USABLE", "USABLE",
            "USABLE_GEOMETRY_OPEN_OBLIGATION",
            "USABLE_GEOMETRY_OPEN_OBLIGATION",
        ])

    def test_P4_and_P5_never_report_their_open_obligations_as_closed(self):
        ctx = _ctx()
        for stage in (REPORT._p4_stage(ctx), REPORT._p5_stage(ctx)):
            self.assertIsNotNone(stage["open_obligation"])
            self.assertNotEqual(stage["verdict"], "USABLE")

    def test_a_microscopic_P1_handoff_is_rejected(self):
        ctx = _ctx()
        ctx["P1_normal_gauged_cayley_norm_upper"] = 2.9e-32
        stage = REPORT._p1_stage(ctx)
        self.assertEqual(stage["verdict"], "NOT_USABLE")
        failed = [c["check"] for c in stage["usability_checks"] if not c["pass"]]
        self.assertIn("handoff family is not microscopic", failed)

    def test_a_shrunk_handoff_domain_is_rejected(self):
        for key, value in (("position_error_norm_upper_m", 1.0),
                           ("velocity_error_norm_upper_mps", 0.1),
                           ("latent_acceleration_error_norm_upper_mps2", 0.1)):
            ctx = _ctx()
            ctx["startup"]["operating_domain"]["physical_handoff_coordinate_bounds"][key] = value
            self.assertEqual(REPORT._p1_stage(ctx)["verdict"], "NOT_USABLE", key)

    def test_a_narrowed_P4_sector_is_rejected(self):
        ctx = _ctx()
        ctx["P4_sector"]["design_full_attitude_angle_rad"] = 0.5
        self.assertEqual(REPORT._p4_stage(ctx)["verdict"], "NOT_USABLE")

    def test_a_cayley_chart_at_or_above_one_is_rejected(self):
        ctx = _ctx()
        ctx["P4_sector"]["design_cayley_norm_upper"] = 1.25
        self.assertEqual(REPORT._p4_stage(ctx)["verdict"], "NOT_USABLE")

    def test_a_shrunk_P5_entrance_is_rejected(self):
        ctx = _ctx()
        ctx["P5_outer"]["declared_P5_entrance"]["gauged_full_attitude_angle_upper_deg"] = 5.0
        self.assertEqual(REPORT._p5_stage(ctx)["verdict"], "NOT_USABLE")

    def test_a_nonzero_outer_word_count_is_rejected(self):
        ctx = _ctx()
        ctx["P5_outer"]["N_outer_words"] = 3
        self.assertEqual(REPORT._p5_stage(ctx)["verdict"], "NOT_USABLE")

    def test_P3_margin_is_labelled_as_a_relative_constant(self):
        stage = REPORT._p3_stage(_ctx())
        self.assertIn("not a nonlinear state radius", stage["interpretation"])

    def test_P4_headline_carries_the_first_accelerometer_budget_gap(self):
        head = REPORT._p4_stage(_ctx())["headline"]["first_accelerometer_sector_budget"]
        self.assertFalse(head["shrinking_the_candidate_angle_alone_can_close_the_budget"])
        self.assertGreater(head["smallest_probe_nuisance_over_budget_ratio"], 1.0)


class ReportValidationTests(unittest.TestCase):
    def _report(self):
        ctx = _ctx()
        stages = [REPORT._p1_stage(ctx), REPORT._p2_stage(ctx), REPORT._p3_stage(ctx),
                  REPORT._p4_stage(ctx), REPORT._p5_stage(ctx)]
        return {
            "schema": REPORT.SCHEMA,
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "filter_changed": False,
            "upstream_pass_fields_trusted_as_usability": False,
            "stages": stages,
            "stages_not_usable": [],
            "open_theorem_obligations": {
                s["stage"]: s["open_obligation"] for s in stages if s["open_obligation"]},
            "all_stage_geometries_usable": True,
            "P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED": False,
            "failures": [],
        }

    def test_a_clean_report_validates(self):
        self.assertEqual(REPORT.validate(self._report()), [])

    def test_complete_proof_cannot_be_claimed_with_open_obligations(self):
        d = self._report()
        d["P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED"] = True
        self.assertIn("complete proof claimed with open obligations", REPORT.validate(d))

    def test_a_failed_usability_check_surfaces_in_validation(self):
        d = deepcopy(self._report())
        d["stages"][0]["usability_checks"][1]["pass"] = False
        self.assertTrue(any("usability check failed" in x for x in REPORT.validate(d)))

    def test_out_of_order_stages_are_rejected(self):
        d = deepcopy(self._report())
        d["stages"] = list(reversed(d["stages"]))
        self.assertIn("stage list is not P1..P5", REPORT.validate(d))


if __name__ == "__main__":
    unittest.main()
