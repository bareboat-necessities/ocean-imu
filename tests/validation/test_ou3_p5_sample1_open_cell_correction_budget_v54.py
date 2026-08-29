import copy
import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51
import ou3_p5_sample1_open_cell_correction_budget_v54 as V54


def _cell(label, corr, admissible):
    spec = V54.V53_OPEN_CELLS[label]
    return {
        "cell": list(spec["cell"]),
        "V53_record_source": spec["record_source"],
        "V53_recorded_correction_radial_upper_rad": spec[
            "correction_radial_upper_rad"],
        "reconstruction_over_V53_recorded_ratio": (
            corr / spec["correction_radial_upper_rad"]),
        "sample1_residual_norm_upper_mps2": 18.0,
        "sample1_combined_source_x_residual_upper_mps2": 1.05,
        "Ktheta_perpendicular_block_upper": 0.76,
        "Ktheta_parallel_block_upper": 0.098,
        "combined_directional_correction_norm_upper_rad": corr,
        "geodesic_admissible_correction_rad": admissible,
        "correction_gap_to_geodesic_target_rad": corr - admissible,
        "correction_fractional_gap": (corr - admissible) / corr,
        "already_inside_geodesic_target": corr <= admissible,
        "required_reductions": {},
    }


def _ok_artifact():
    cells = {
        "first_open": _cell("first_open", 1.9595039241752017, 1.9389861206167838),
        "worst": _cell("worst", 3.1985263367953887, 1.257263970049669),
        "retired_witness": _cell(
            "retired_witness", V54.V51_RETIRED_WITNESS_CORRECTION_RAD,
            2.0308497552113565),
    }
    return {
        "schema": V54.SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54",
        "source_generated_not_trajectory_fit": True,
        "exact_monotone_corner_enclosure_used": True,
        "required_reductions_are_diagnostics_not_enclosures": True,
        "geodesic_branch_only": True,
        "source_replay_used": False,
        "filter_changed": False,
        "cell_reported_closed_or_unreachable_here": False,
        "deployed_correction_limit_increased": False,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "deployed_correction_limit_rad": 6.0,
        "q_target": V54.Q_TARGET,
        "open_cells": cells,
        "nearest_open_cell": "first_open",
        "nearest_open_cell_fractional_gap": 0.010470917309876947,
        "P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54": "PASS",
        "next_obligation": "REFINE_NOMINAL_SAMPLE1_RESIDUAL_AT_THE_NEAREST_OPEN_Q8_CELL",
        "failures": [],
    }


class Sample1OpenCellCorrectionBudgetV54Tests(unittest.TestCase):
    def test_targets_and_recorded_cells_are_frozen(self):
        self.assertEqual(V54.SCHEMA, 5400)
        self.assertEqual(V54.Q_TARGET, 8.0)
        self.assertEqual(V54.V53_OPEN_CELLS["first_open"]["cell"], (0, 2, 23))
        self.assertEqual(V54.V53_OPEN_CELLS["worst"]["cell"], (23, 18, 4))
        self.assertEqual(V54.V53_OPEN_CELLS["retired_witness"]["cell"], (0, 0, 23))
        self.assertEqual(V54.V51_RETIRED_WITNESS_CORRECTION_RAD,
                         1.7313776836494923)
        self.assertEqual(
            V54.V53_OPEN_CELLS["retired_witness"][
                "sample1_current_cayley_norm_upper"], V51.V41_Q_CURRENT)

    def test_admissible_correction_matches_the_triangle(self):
        # phi_c + phi_d = 2 atan(Q_TARGET/2) at the admissible correction.
        for q in (0.0, 0.6415212986499801, 1.674977608400414, 3.0):
            adm = V54._admissible_correction(q)
            phi_c = 2.0 * math.atan(0.5 * q)
            self.assertAlmostEqual(phi_c + adm, 2.0 * math.atan(4.0), places=12)
        self.assertGreater(V54._admissible_correction(0.6415212986499801),
                           V54._admissible_correction(1.674977608400414))

    def test_admissible_correction_fails_closed_on_bad_input(self):
        with self.assertRaises(ValueError):
            V54._admissible_correction(-1.0)
        with self.assertRaises(ValueError):
            V54._admissible_correction(float("nan"))

    def test_required_reduction_reaches_the_target_when_applied(self):
        k_perp, k_par, rho, rho_x = 0.761921, 0.097729, 18.3283, 1.05144
        target = 1.9389861206167838
        out = V54._required_reductions(k_perp=k_perp, k_par=k_par, rho=rho,
                                       rho_x=rho_x, target=target)
        for name in ("rho", "rho_x", "k_parallel", "k_perpendicular"):
            row = out[name]
            self.assertTrue(row["reachable_alone"], name)
            kp, kq, rr, rx = k_perp, k_par, rho, rho_x
            if name == "rho":
                rr = row["required"]
            elif name == "rho_x":
                rx = row["required"]
            elif name == "k_parallel":
                kq = row["required"]
            else:
                kp = row["required"]
            got = math.sqrt(kp * kp * rx * rx + kq * kq * (rr * rr - rx * rx))
            self.assertAlmostEqual(got, target, places=9, msg=name)

    def test_unreachable_factor_is_reported_not_invented(self):
        # k_par * rho alone already exceeds the target, so zeroing the
        # perpendicular factors cannot reach it.
        out = V54._required_reductions(k_perp=0.6038, k_par=0.16265,
                                       rho=19.3024, rho_x=1.05144,
                                       target=1.257263970049669)
        self.assertFalse(out["k_perpendicular"]["reachable_alone"])
        self.assertFalse(out["rho_x"]["reachable_alone"])
        self.assertTrue(out["rho"]["reachable_alone"])
        self.assertTrue(out["k_parallel"]["reachable_alone"])

    def test_required_reductions_fail_closed_on_bad_input(self):
        with self.assertRaises(ValueError):
            V54._required_reductions(k_perp=-1.0, k_par=0.1, rho=1.0,
                                     rho_x=0.5, target=0.5)
        with self.assertRaises(ValueError):
            V54._required_reductions(k_perp=1.0, k_par=0.1, rho=0.5,
                                     rho_x=1.0, target=0.5)

    def test_validate_accepts_a_well_formed_artifact(self):
        self.assertEqual(V54.validate(_ok_artifact()), [])

    def test_validate_rejects_a_drifted_or_promoted_artifact(self):
        drifted = copy.deepcopy(_ok_artifact())
        drifted["open_cells"]["retired_witness"][
            "combined_directional_correction_norm_upper_rad"] = 1.5
        self.assertTrue(V54.validate(drifted))

        moved = copy.deepcopy(_ok_artifact())
        moved["open_cells"]["worst"]["cell"] = [1, 2, 3]
        self.assertTrue(V54.validate(moved))

        dropped = copy.deepcopy(_ok_artifact())
        del dropped["open_cells"]["worst"]
        self.assertTrue(V54.validate(dropped))

        inconsistent = copy.deepcopy(_ok_artifact())
        inconsistent["open_cells"]["first_open"][
            "sample1_combined_source_x_residual_upper_mps2"] = 99.0
        self.assertTrue(V54.validate(inconsistent))

        overclaimed = copy.deepcopy(_ok_artifact())
        overclaimed["cell_reported_closed_or_unreachable_here"] = True
        self.assertTrue(V54.validate(overclaimed))

        not_diagnostic = copy.deepcopy(_ok_artifact())
        not_diagnostic["required_reductions_are_diagnostics_not_enclosures"] = False
        self.assertTrue(V54.validate(not_diagnostic))

        promoted = copy.deepcopy(_ok_artifact())
        promoted["N_H_words_set_here"] = True
        self.assertTrue(V54.validate(promoted))

        limit = copy.deepcopy(_ok_artifact())
        limit["deployed_correction_limit_rad"] = 9.0
        self.assertTrue(V54.validate(limit))


if __name__ == "__main__":
    unittest.main()
