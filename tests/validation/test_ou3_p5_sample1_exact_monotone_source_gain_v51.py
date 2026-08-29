import math
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51

G = 9.80665


def _box(lo, hi):
    return Interval.outward_bounds(lo, hi)


def _exact_sample(t, p, r, g):
    """Reference values of the eight block rationals at one parameter point."""
    D = g * g * t + p + r
    return {
        "a": t * (p + r) / D,
        "c0": -(g * t * p / D),
        "b": p * (g * g * t + r) / D,
        "bz": p * r / (p + r),
        "det_first": t * p * r / D,
        "ktheta": g * t / D,
        "kaw_t": p / D,
        "kz": p / (p + r),
    }


class Sample1ExactMonotoneSourceGainV51Tests(unittest.TestCase):
    def test_authoritative_witness_and_targets_are_frozen(self):
        self.assertEqual(V51.WITNESS, (0, 0, 23))
        self.assertEqual(V51.SCHEMA, 5100)
        self.assertEqual(V51.Q_TARGET, 8.0)
        self.assertEqual(V51.V41_Q_CURRENT, 0.6415212986499801)
        self.assertEqual(V51.V41_Q_POST, 8.344528951460543)

    def test_exact_block_is_inside_the_parent_block(self):
        t = _box(1.0e-3, 2.0e-3)
        p = _box(2.0e-3, 4.0e-3)
        r = _box(8.0e-4, 9.0e-4)
        parent = V51._first_block(t=t, p=p, r=r, g=G, exact=False)
        exact = V51._first_block(t=t, p=p, r=r, g=G, exact=True)
        for key in parent:
            self.assertGreaterEqual(exact[key].lo, parent[key].lo, key)
            self.assertLessEqual(exact[key].hi, parent[key].hi, key)

    def test_exact_block_encloses_every_sampled_parameter_point(self):
        rng = random.Random(20260829)
        for _ in range(64):
            tl = rng.uniform(1.0e-4, 5.0e-3)
            pl = rng.uniform(1.0e-4, 5.0e-2)
            rl = rng.uniform(1.0e-5, 5.0e-3)
            t = _box(tl, tl * rng.uniform(1.0, 4.0))
            p = _box(pl, pl * rng.uniform(1.0, 4.0))
            r = _box(rl, rl * rng.uniform(1.0, 4.0))
            exact = V51._first_block(t=t, p=p, r=r, g=G, exact=True)
            for _s in range(12):
                sample = _exact_sample(
                    rng.uniform(t.lo, t.hi), rng.uniform(p.lo, p.hi),
                    rng.uniform(r.lo, r.hi), G)
                for key, value in sample.items():
                    self.assertGreaterEqual(
                        value, exact[key].lo,
                        f"{key} below refined lower bound")
                    self.assertLessEqual(
                        value, exact[key].hi,
                        f"{key} above refined upper bound")

    def test_axial_aw_gain_upper_drops_below_one(self):
        # p/(p+r) < 1 identically, but the parent interval expression can
        # exceed one.  The corner enclosure must not.
        t = _box(1.225e-3, 1.225e-3)
        p = _box(2.4261138354980992e-3, 3.4696585695800247e-3)
        r = _box(8.666185998411353e-4, 8.666185998411355e-4)
        parent = V51._first_block(t=t, p=p, r=r, g=G, exact=False)
        exact = V51._first_block(t=t, p=p, r=r, g=G, exact=True)
        self.assertGreater(parent["kz"].hi, 1.0)
        self.assertLess(exact["kz"].hi, 1.0)
        self.assertLess(exact["kz"].hi, parent["kz"].hi)

    def test_monotone_intersection_fails_closed_on_a_disjoint_parent(self):
        good = V51._monotone(Interval(1.0, 1.0), Interval(2.0, 2.0),
                             Interval(0.0, 3.0), name="ok")
        self.assertEqual((good.lo, good.hi), (1.0, 2.0))
        with self.assertRaises(RuntimeError):
            V51._monotone(Interval(5.0, 5.0), Interval(6.0, 6.0),
                          Interval(0.0, 1.0), name="bad")

    def test_geodesic_composition_rejects_out_of_range_corrections(self):
        self.assertEqual(V51._geodesic_q_upper(-1.0), math.inf)
        self.assertEqual(V51._geodesic_q_upper(4.0), math.inf)
        self.assertEqual(V51._geodesic_q_upper(float("nan")), math.inf)
        parent_q = V51._geodesic_q_upper(2.0466720610769817)
        refined_q = V51._geodesic_q_upper(1.7313776836494923)
        self.assertGreater(parent_q, V51.Q_TARGET)
        self.assertLess(refined_q, V51.Q_TARGET)
        self.assertLess(refined_q, parent_q)

    def test_provenance_rejects_a_drifted_parent_reconstruction(self):
        ok = dict(V51.PARENT_WITNESS)
        self.assertEqual(V51._provenance(ok), [])
        drifted = dict(ok)
        drifted["sample1_force_norm_upper_mps2"] = 21.0
        self.assertEqual(len(V51._provenance(drifted)), 1)

        shifted = dict(ok)
        shifted["first_axial_residual_mps2"] = [10.219984799344509, 11.15]
        self.assertEqual(len(V51._provenance(shifted)), 1)

        self.assertEqual(len(V51._provenance({})), len(V51.PARENT_WITNESS))

    def test_validate_rejects_a_promoted_or_widened_artifact(self):
        base = {
            "schema": V51.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51",
            "source_generated_not_trajectory_fit": True,
            "exact_monotone_corner_enclosure_used": True,
            "parent_enclosure_retained_as_intersection": True,
            "archived_parent_witness_reproduced": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "full_source_cell0_cover_lifted_here": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_established_here": False,
            "V41_first_survivor_row": list(V51.WITNESS),
            "deployed_correction_limit_rad": 6.0,
            "q_target": V51.Q_TARGET,
            "parent_witness_chain": {
                "first_block": {"kz": [0.55, 1.06]},
                "combined_directional_correction_norm_upper_rad": 2.0466720610769817,
            },
            "exact_monotone_witness_chain": {
                "first_block": {"kz": [0.73, 0.81]},
                "combined_directional_correction_norm_upper_rad": 1.7313776836494923,
            },
            "witness_comparison": {
                "parent_geodesic_q_upper": 8.277775079246437,
                "refined_geodesic_q_upper": 4.8010333986449245,
                "authoritative_witness_closed_by_refined_correction": True,
            },
            "authoritative_witness_closed": True,
            "P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51": "PASS",
            "next_obligation": "LIFT_EXACT_MONOTONE_SOURCE_GAIN_OVER_FULL_V41_SOURCE_CELL0_COVER",
            "failures": [],
        }
        self.assertEqual(V51.validate(base), [])

        promoted = dict(base, N_H_words_set_here=True)
        self.assertTrue(V51.validate(promoted))

        lifted = dict(base, full_source_cell0_cover_lifted_here=True)
        self.assertTrue(V51.validate(lifted))

        escaped = dict(base)
        escaped["exact_monotone_witness_chain"] = dict(
            base["exact_monotone_witness_chain"],
            first_block={"kz": [0.40, 0.81]})
        self.assertTrue(V51.validate(escaped))

        widened = dict(base)
        widened["exact_monotone_witness_chain"] = dict(
            base["exact_monotone_witness_chain"],
            combined_directional_correction_norm_upper_rad=3.0)
        self.assertTrue(V51.validate(widened))

        inconsistent = dict(base)
        inconsistent["witness_comparison"] = dict(
            base["witness_comparison"],
            authoritative_witness_closed_by_refined_correction=False)
        self.assertTrue(V51.validate(inconsistent))

        unreproduced = dict(base, archived_parent_witness_reproduced=False)
        self.assertTrue(V51.validate(unreproduced))

        limit = dict(base, deployed_correction_limit_rad=9.0)
        self.assertTrue(V51.validate(limit))


if __name__ == "__main__":
    unittest.main()
