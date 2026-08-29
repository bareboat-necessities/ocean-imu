import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_zero_perturbation_barrier_v50 as V50


def _terms(**over):
    args = dict(h_row_norm=20.0, h_norm=21.0, hp_row_norm=1.0e-2,
                hp_row_theta_norm=9.0e-3, p_norm=8.0e-3, p_theta_norm=7.0e-3,
                c_theta_norm=1.5e-2, dP=3.0e-10, dH=1.0e-7)
    args.update(over)
    return V50._row_terms(**args)


class Sample1ZeroPerturbationBarrierV50Tests(unittest.TestCase):
    def test_authoritative_witness_and_targets_are_frozen(self):
        self.assertEqual(V50.WITNESS, (0, 0, 23))
        self.assertEqual(V50.SCHEMA, 5000)
        self.assertEqual(V50.Q_TARGET, 8.0)

    def test_attitude_supported_terms_never_exceed_their_parent(self):
        d = _terms()
        self.assertTrue(d["attitude_supported_terms_never_exceed_parent"])
        for key in d["term_order"]:
            self.assertLessEqual(d["refined_terms"][key],
                                 d["V34_parent_terms"][key])
        self.assertLessEqual(d["attitude_supported_row_candidate_upper"],
                             d["V34_parent_row_candidate_upper"])

    def test_only_the_three_delta_h_terms_are_restricted(self):
        d = _terms()
        unchanged = ("hi_dP_H", "dH_dP_H", "hi_dP_dH", "dH_dP_dH")
        for key in unchanged:
            self.assertEqual(d["attitude_supported_terms"][key],
                             d["V34_parent_terms"][key])
        for key in ("dH_P_H", "hiP_dH", "dH_P_dH"):
            self.assertLess(d["attitude_supported_terms"][key],
                            d["V34_parent_terms"][key])

    def test_a_worse_restricted_factor_cannot_raise_a_term(self):
        # A hypothetical attitude-restricted factor larger than its parent must
        # be discarded by the per-term minimum rather than widen the bound.
        d = _terms(c_theta_norm=1.0e6, hp_row_theta_norm=1.0e6,
                   p_theta_norm=1.0e6)
        self.assertTrue(d["attitude_supported_terms_never_exceed_parent"])
        self.assertEqual(d["attitude_supported_row_candidate_upper"],
                         d["V34_parent_row_candidate_upper"])

    def test_row_terms_fail_closed_on_invalid_input(self):
        with self.assertRaises(ValueError):
            _terms(dP=-1.0)
        with self.assertRaises(ValueError):
            _terms(h_norm=math.inf)

    def test_reduced_covariance_decomposition_reproduces_its_parent(self):
        FULL = V50.FULL
        dPplus, dirterm, eps, s_part, dhi = (
            4.2e-11, 8.15e-11, 1.48e-10, 5.1e-12, 0.059)
        tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
        psd = FULL.up(FULL.up(FULL.up(tnom * tnom * dPplus) + dirterm) + eps)
        vr = {
            "first_posterior_covariance_perturbation_upper": dPplus,
            "reset_gauge_transform_perturbation_upper": dirterm,
            "sample1_reduced_covariance_PSD_perturbation_upper": psd,
            "S_reduced_covariance_perturbation_upper": s_part,
            "total_reduced_covariance_perturbation_upper": FULL.up(psd + s_part),
        }
        d = V50._reduced_covariance_terms(vr=vr, dhi=dhi, eps=eps)
        self.assertTrue(d["PSD_decomposition_reproduces_certified_parent"])
        self.assertTrue(d["total_decomposition_reproduces_certified_parent"])
        self.assertIn(d["dominant_term"], V50._NEXT_BY_DP_TERM)
        self.assertEqual(d["dominant_term"],
                         "sample1_prediction_attitude_epsilon_upper")

    def test_reduced_covariance_decomposition_detects_a_broken_parent(self):
        vr = {
            "first_posterior_covariance_perturbation_upper": 1.0e-11,
            "reset_gauge_transform_perturbation_upper": 1.0e-11,
            "sample1_reduced_covariance_PSD_perturbation_upper": 1.0e-11,
            "S_reduced_covariance_perturbation_upper": 1.0e-12,
            "total_reduced_covariance_perturbation_upper": 1.0e-11,
        }
        d = V50._reduced_covariance_terms(vr=vr, dhi=0.05, eps=1.0e-10)
        self.assertFalse(d["total_decomposition_reproduces_certified_parent"])

    def test_zeroed_caps_zero_every_perturbation_radius(self):
        caps = {k: 1.0 for k in V50._ZEROED_CAP_KEYS}
        caps["V31_parent_yz_norm_upper_rad"] = 3.0
        out = V50._zeroed_caps(caps)
        for key in V50._ZEROED_CAP_KEYS:
            self.assertEqual(out[key], 0.0)
        self.assertEqual(out["V31_parent_yz_norm_upper_rad"], 3.0)
        self.assertIsNot(out, caps)
        self.assertEqual(caps[V50._ZEROED_CAP_KEYS[0]], 1.0)

    def test_zeroed_caps_fail_closed_on_a_missing_radius(self):
        caps = {k: 1.0 for k in V50._ZEROED_CAP_KEYS[:-1]}
        with self.assertRaises(KeyError):
            V50._zeroed_caps(caps)

    def test_angle_diagnostics_measure_the_geodesic_shortfall(self):
        d = V50._angle_diagnostics(q_current=0.6415212986499801,
                                   radial_lower=1.1782765791511944,
                                   radial_upper=2.050326092645528,
                                   cap_total=3.654031568545583e-3)
        self.assertGreater(d["principal_angle_reduction_needed_rad"], 0.0)
        self.assertLess(d["cap_share_of_needed_reduction"], 1.0)
        self.assertAlmostEqual(d["q_target_principal_angle_rad"],
                               2.0 * math.atan(4.0))
        # Removing the whole certified perturbation budget from the radial
        # bound still leaves the geodesic branch above the q<8 target.
        self.assertLess(d["cap_free_correction_principal_angle_rad"],
                        d["correction_principal_angle_rad"])
        self.assertTrue(d["cap_free_geodesic_q_still_at_or_above_target"])
        self.assertGreater(d["cap_free_geodesic_q_diagnostic_upper"],
                           V50.Q_TARGET)

    def test_angle_diagnostics_fail_closed_on_invalid_input(self):
        with self.assertRaises(ValueError):
            V50._angle_diagnostics(q_current=-1.0, radial_lower=1.0,
                                   radial_upper=2.0, cap_total=0.0)

    def test_validate_rejects_a_promoted_or_unbound_artifact(self):
        base = {
            "schema": V50.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50",
            "source_generated_not_trajectory_fit": True,
            "attitude_supported_Jacobian_perturbation_used": True,
            "V12D_full_DeltaC_parent_retained": True,
            "V12D_full_DeltaS_parent_retained_as_intersection": True,
            "V34_seven_term_row_expansion_retained": True,
            "temporary_V48_hooks_restored": True,
            "zero_perturbation_run_is_route_audit_not_filter_claim": True,
            "source_replay_used": False,
            "filter_changed": False,
            "failed_V33_row_candidate_promoted": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_established_here": False,
            "V41_first_survivor_row": list(V50.WITNESS),
            "deployed_correction_limit_rad": 6.0,
            "q_target": V50.Q_TARGET,
            "attitude_supported_row_detail": {
                "V34_first_measurement_row_DeltaS_intersected_upper": 2.0e-7,
                "attitude_supported_row_DeltaS_intersected_upper": 1.0e-7,
                "DeltaC_parent_is_exact_attitude_supported_expansion": True,
            },
            "reduced_covariance_perturbation_decomposition": {
                "total_decomposition_reproduces_certified_parent": True,
            },
            "zero_perturbation_barrier": {
                "zero_perturbation_best_q_upper": 8.29,
                "zero_perturbation_joint_box_nonempty": True,
                "barrier_established": True,
                "perturbation_route_can_close_authoritative_witness": False,
            },
            "P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50": "PASS",
            "next_obligation": V50._BARRIER_OBLIGATION,
            "failures": [],
        }
        self.assertEqual(V50.validate(base), [])

        promoted = dict(base, N_H_words_set_here=True)
        self.assertTrue(V50.validate(promoted))

        widened = dict(base)
        widened["attitude_supported_row_detail"] = dict(
            base["attitude_supported_row_detail"],
            attitude_supported_row_DeltaS_intersected_upper=3.0e-7)
        self.assertTrue(V50.validate(widened))

        unbound = dict(base)
        unbound["reduced_covariance_perturbation_decomposition"] = {
            "total_decomposition_reproduces_certified_parent": False}
        self.assertTrue(V50.validate(unbound))

        empty_box = dict(base)
        empty_box["zero_perturbation_barrier"] = dict(
            base["zero_perturbation_barrier"],
            zero_perturbation_joint_box_nonempty=False)
        self.assertTrue(V50.validate(empty_box))

        inconsistent = dict(base)
        inconsistent["zero_perturbation_barrier"] = dict(
            base["zero_perturbation_barrier"],
            perturbation_route_can_close_authoritative_witness=True)
        self.assertTrue(V50.validate(inconsistent))

        limit = dict(base, deployed_correction_limit_rad=9.0)
        self.assertTrue(V50.validate(limit))

        target = dict(base, q_target=9.0)
        self.assertTrue(V50.validate(target))


if __name__ == "__main__":
    unittest.main()
