from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_domain_retention as D  # noqa: E402


class DomainRetentionTests(unittest.TestCase):
    def test_declared_radii_come_from_the_operating_domain_file(self):
        """Radii are read from the declared domain, never hard-coded here."""
        radii = D.declared_radii()
        self.assertEqual(set(radii), {name for name, _, _ in D.GROUPS})
        self.assertAlmostEqual(radii["attitude"], 2*np.tan(np.deg2rad(15.)), places=12)
        self.assertEqual(radii["accelerometer_bias"], .4)
        self.assertEqual(radii["velocity"], 5.)
        self.assertEqual(radii["integral_displacement"], 300.)
        for value in radii.values():
            self.assertGreater(value, 0.)

    def test_groups_tile_the_twenty_one_error_coordinates_in_order(self):
        """The seven groups partition the 21 error coordinates without a gap."""
        self.assertEqual([offset for _, offset, _ in D.GROUPS], list(range(0, 21, 3)))

    def example(self):
        """A small random word shared by the product-set tests."""
        rng = np.random.default_rng(2026)
        blocks = {name: rng.normal(size=(3, 3))/3 for name, _, _ in D.GROUPS}
        radii = D.declared_radii()
        return blocks, radii, rng.normal(size=3)/10

    def test_attained_bound_never_exceeds_its_certified_enclosure(self):
        """The attained functional is a lower bound on the same supremum."""
        blocks, radii, forcing = self.example()
        lower = D.attained_lower_bound(blocks, radii, forcing)
        upper = D.certified_upper_bound(blocks, radii, forcing)
        self.assertLessEqual(lower, upper*(1+1e-12))
        self.assertGreater(lower, 0.)

    def test_attained_bound_dominates_every_sampled_admissible_state(self):
        """No sampled admissible initial state reaches past the attained bound."""
        blocks, radii, forcing = self.example()
        lower = D.attained_lower_bound(blocks, radii, forcing)
        rng = np.random.default_rng(7)
        for _ in range(500):
            reached = forcing*rng.uniform(-1., 1.)
            for name, block in blocks.items():
                x = rng.normal(size=3)
                reached = reached + radii[name]*block@(x/np.linalg.norm(x))
            self.assertLessEqual(np.linalg.norm(reached), lower*(1+1e-9))

    def test_source_contributions_reproduce_the_certified_bound(self):
        """The per-source shares sum back to the subadditive bound they explain."""
        blocks, radii, forcing = self.example()
        shares = D.source_contributions(blocks, radii, forcing)
        self.assertEqual(shares["forcing_template"], float(np.linalg.norm(forcing)))
        self.assertAlmostEqual(sum(shares.values()),
                               D.certified_upper_bound(blocks, radii, forcing), places=12)

    def test_retention_is_positively_homogeneous_in_the_declared_radii(self):
        """Scaling every radius and the template together leaves the ratios fixed."""
        rng = np.random.default_rng(11)
        transitions = [rng.normal(size=(21, 21))/4 for _ in range(3)]
        responses = [rng.normal(size=21)/8 for _ in range(3)]
        radii = D.declared_radii()
        base = D.retention(transitions, responses, radii)
        scaled = D.retention(transitions, [2*r for r in responses],
                             {k: 2*v for k, v in radii.items()})
        for name, item in base.items():
            self.assertAlmostEqual(scaled[name]["certified_retention_ratio"],
                                   item["certified_retention_ratio"], places=10)
            self.assertEqual(scaled[name]["certified_prefix_index"],
                             item["certified_prefix_index"])

    def test_attained_bound_is_taken_over_every_prefix(self):
        """The loosest subadditive prefix need not be the one reached furthest.

        The first prefix sends the seven declared balls onto three mutually
        orthogonal output axes, so no single direction collects them all and
        its subadditive sum 7 overstates the reached 17**.5. The second sends
        one ball along one axis, where the two bounds agree at 5. The
        certified maximum is therefore the first prefix and the attained
        maximum the second, and selecting on the certified bound alone would
        report 17**.5 and miss the larger 5.
        """
        radii = D.declared_radii()
        spread, concentrated = np.zeros((21, 21)), np.zeros((21, 21))
        for position, (name, start, _) in enumerate(D.GROUPS):
            axis = np.zeros(3)
            axis[position % 3] = 1.
            spread[0:3, start:start+3] = np.outer(axis, [1., 0., 0.])/radii[name]
        concentrated[0:3, 12:15] = 5*np.outer([1., 0., 0.], [1., 0., 0.])/radii["integral_displacement"]
        transitions = [spread, concentrated]
        responses = [np.zeros(21), np.zeros(21)]
        result = D.retention(transitions, responses, radii)["attitude"]
        self.assertEqual(result["certified_prefix_index"], 0)
        self.assertEqual(result["attained_prefix_index"], 1)
        self.assertAlmostEqual(result["certified_upper_bound"], 7., places=9)
        self.assertAlmostEqual(result["attained_lower_bound"], 5., places=9)

    def test_a_boundary_ratio_is_not_reported_as_a_definite_violation(self):
        """A group exactly at its radius sits inside the shared margin."""
        radii = D.declared_radii()
        blocks = {"attitude": np.eye(3)}
        active = {"attitude": radii["attitude"]}
        forcing = np.zeros(3)
        attained = D.attained_lower_bound(blocks, active, forcing)
        self.assertLessEqual(attained/radii["attitude"], 1.+D.MARGIN)
        self.assertGreater(D.MARGIN, 0.)

    def test_restricting_the_initial_set_can_only_lower_the_bound(self):
        """Dropping initial balls can never raise the reachable excursion."""
        rng = np.random.default_rng(13)
        transitions = [rng.normal(size=(21, 21))/4 for _ in range(3)]
        responses = [rng.normal(size=21)/8 for _ in range(3)]
        radii = D.declared_radii()
        full = D.retention(transitions, responses, radii)
        bias = D.retention(transitions, responses, radii, keep={"accelerometer_bias"})
        for name, item in bias.items():
            self.assertLessEqual(item["certified_upper_bound"],
                                 full[name]["certified_upper_bound"]*(1+1e-12))


    def ellipsoid_example(self):
        """A small random word and covariance shared by the ellipsoid tests."""
        rng = np.random.default_rng(31)
        factor = rng.normal(size=(21, 21))/8
        covariance = factor@factor.T + np.eye(21)*1e-3
        transitions = [np.eye(21), rng.normal(size=(21, 21))/10]
        responses = [np.zeros(21), rng.normal(size=21)/100]
        return transitions, responses, covariance

    def test_ellipsoid_excursion_dominates_every_sampled_admissible_state(self):
        """The ellipsoid image needs no subadditive step, so it is exact."""
        transitions, responses, covariance = self.ellipsoid_example()
        radii = D.declared_radii()
        out = D.ellipsoid_retention(transitions, responses, radii, covariance)
        root = np.linalg.cholesky(covariance)
        rng = np.random.default_rng(37)
        for name, offset, _ in D.GROUPS:
            rows = slice(offset, offset+3)
            reach = out["groups"][name]["excursion_at_one_sigma"]
            index = out["groups"][name]["one_sigma_prefix_index"]
            for _ in range(300):
                y = rng.normal(size=21)
                x = root @ (y/np.linalg.norm(y))
                state = transitions[index] @ x + responses[index]*rng.uniform(-1., 1.)
                self.assertLessEqual(np.linalg.norm(state[rows]), reach*(1+1e-9))

    def test_critical_level_scales_inversely_with_the_initial_covariance(self):
        """Quadrupling the covariance halves every critical level."""
        transitions, responses, covariance = self.ellipsoid_example()
        radii = D.declared_radii()
        base = D.ellipsoid_retention(transitions, responses, radii, covariance)
        scaled = D.ellipsoid_retention(transitions, responses, radii, 4*covariance)
        for name, item in base["groups"].items():
            self.assertAlmostEqual(scaled["groups"][name]["critical_initial_sigma_level"],
                                   item["critical_initial_sigma_level"]/2., places=9)

    def test_a_group_the_forcing_alone_evicts_admits_no_level(self):
        """Forcing past a declared radius leaves no admissible level at all."""
        transitions, responses, covariance = self.ellipsoid_example()
        radii = D.declared_radii()
        evicted = [np.zeros(21), np.zeros(21)]
        evicted[1][3:6] = 2*radii["gyro_bias"]
        out = D.ellipsoid_retention(transitions, evicted, radii, covariance)
        self.assertEqual(out["groups"]["gyro_bias"]["critical_initial_sigma_level"], 0.)
        self.assertEqual(out["limiting_group"], "gyro_bias")
        self.assertFalse(out["P4_PASS"])

    def test_a_zero_gain_prefix_never_limits_the_level(self):
        """A prefix with no gain places no bound on the level."""
        radii = D.declared_radii()
        covariance = np.eye(21)
        transitions = [np.zeros((21, 21)), np.eye(21)]
        responses = [np.zeros(21), np.zeros(21)]
        out = D.ellipsoid_retention(transitions, responses, radii, covariance)
        for name, item in out["groups"].items():
            self.assertAlmostEqual(item["critical_initial_sigma_level"],
                                   radii[name], places=9)


    def test_optimal_storage_attains_the_squared_spectral_radius(self):
        """The eigenbasis witness meets the infimum over all metrics."""
        rng = np.random.default_rng(101)
        for _ in range(5):
            transition = rng.normal(size=(6, 6))/3
            out = D.optimal_quadratic_storage(transition)
            spectral = max(abs(np.linalg.eigvals(transition)))
            self.assertAlmostEqual(out["spectral_radius"], spectral, places=10)
            self.assertAlmostEqual(out["best_achievable_complete_word_ratio"],
                                   spectral**2, places=10)
            self.assertAlmostEqual(out["witness_metric_achieves"],
                                   spectral**2, places=6)

    def test_no_metric_beats_the_spectral_radius(self):
        """Random metrics never fall below the reported infimum."""
        rng = np.random.default_rng(103)
        transition = rng.normal(size=(6, 6))/3
        floor = D.optimal_quadratic_storage(
            transition)["best_achievable_complete_word_ratio"]
        for _ in range(200):
            a = rng.normal(size=(6, 6))
            root = np.linalg.cholesky(a@a.T + np.eye(6))
            ratio = np.linalg.norm(
                root.T @ transition @ np.linalg.inv(root.T), 2)**2
            self.assertGreaterEqual(ratio, floor*(1-1e-9))

    def test_a_unit_eigenvalue_admits_no_contracting_storage(self):
        """An unobserved state puts the spectral radius exactly at one."""
        transition = np.diag([.5, .25, 1.])
        out = D.optimal_quadratic_storage(transition)
        self.assertEqual(out["spectral_radius"], 1.)
        self.assertFalse(out["a_contracting_quadratic_storage_exists"])
        self.assertFalse(out["P4_PASS"])

    def test_similarity_leaves_the_infimum_unchanged(self):
        """The bound is a similarity invariant, not a coordinate artefact."""
        rng = np.random.default_rng(107)
        transition = rng.normal(size=(6, 6))/3
        basis = rng.normal(size=(6, 6)) + 3*np.eye(6)
        moved = np.linalg.solve(basis, transition) @ basis
        self.assertAlmostEqual(
            D.optimal_quadratic_storage(moved)["spectral_radius"],
            D.optimal_quadratic_storage(transition)["spectral_radius"], places=9)


    def cascade_word(self, coupling):
        """A two-step word whose bias couples back into motion by `coupling`."""
        rng = np.random.default_rng(211)
        factor = np.zeros((21, 21))
        factor[:18, :18] = np.linalg.qr(rng.normal(size=(18, 18)))[0]*.9
        factor[:18, 18:] = rng.normal(size=(18, 3))/10
        factor[18:, 18:] = np.eye(3)*.999
        factor[18:, :18] = coupling
        steps = [{"A": factor}, {"A": factor}]
        return steps, factor@factor, np.zeros(21)

    def test_zero_coupling_is_reported_as_an_exact_cascade(self):
        """A vanishing bias<-motion block in every factor is a cascade."""
        steps, transition, response = self.cascade_word(np.zeros((3, 18)))
        out = D.cascade_structure(steps, transition, response)
        self.assertTrue(out["exact_cascade"])
        self.assertEqual(out["per_step_bias_from_motion_defect"], 0.)
        self.assertTrue(out["cascade_iss_applies"])
        self.assertAlmostEqual(out["iss_gain"], 1/(1-out["motion_weighted_norm"]),
                               places=9)
        self.assertFalse(out["P4_PASS"])

    def test_nonzero_coupling_is_not_reported_as_a_cascade(self):
        """Any per-step coupling disqualifies the exact-cascade claim."""
        coupling = np.zeros((3, 18))
        coupling[0, 0] = 1e-6
        steps, transition, response = self.cascade_word(coupling)
        out = D.cascade_structure(steps, transition, response)
        self.assertFalse(out["exact_cascade"])
        self.assertFalse(out["cascade_iss_applies"])
        self.assertGreater(out["per_step_bias_from_motion_defect"], 0.)

    def test_endpoint_cancellation_does_not_pass_as_a_cascade(self):
        """Per-step factors decide, so endpoint cancellation cannot hide."""
        a = np.eye(21)
        a[18, 0] = 1.
        b = np.eye(21)
        b[18, 0] = -1.
        steps = [{"A": a}, {"A": b}]
        product = b @ a
        self.assertEqual(np.linalg.norm(product[18:, :18], 2), 0.)
        out = D.cascade_structure(steps, product, np.zeros(21))
        self.assertFalse(out["exact_cascade"])

    def test_iss_limit_scales_with_the_bias_ball(self):
        """The bounded-bias limit is linear in the ball it is given."""
        steps, transition, response = self.cascade_word(np.zeros((3, 18)))
        one = D.cascade_structure(steps, transition, response, bias_ball=.4)
        two = D.cascade_structure(steps, transition, response, bias_ball=.8)
        gain = one["bias_to_motion_weighted_gain"]/(1-one["motion_weighted_norm"])
        self.assertAlmostEqual(two["iss_limit"]-one["iss_limit"], gain*.4, places=9)


if __name__ == "__main__":
    unittest.main()
