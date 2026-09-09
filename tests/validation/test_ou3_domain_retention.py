from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_domain_retention as D  # noqa: E402


class DomainRetentionTests(unittest.TestCase):
    def test_declared_radii_come_from_the_operating_domain_file(self):
        radii = D.declared_radii()
        self.assertEqual(set(radii), {name for name, _, _ in D.GROUPS})
        self.assertAlmostEqual(radii["attitude"], 2*np.tan(np.deg2rad(15.)), places=12)
        self.assertEqual(radii["accelerometer_bias"], .4)
        self.assertEqual(radii["velocity"], 5.)
        self.assertEqual(radii["integral_displacement"], 300.)
        for value in radii.values():
            self.assertGreater(value, 0.)

    def test_groups_tile_the_twenty_one_error_coordinates_in_order(self):
        self.assertEqual([offset for _, offset, _ in D.GROUPS], list(range(0, 21, 3)))

    def example(self):
        rng = np.random.default_rng(2026)
        blocks = {name: rng.normal(size=(3, 3))/3 for name, _, _ in D.GROUPS}
        radii = D.declared_radii()
        return blocks, radii, rng.normal(size=3)/10

    def test_attained_bound_never_exceeds_its_certified_enclosure(self):
        blocks, radii, forcing = self.example()
        lower = D.attained_lower_bound(blocks, radii, forcing)
        upper = D.certified_upper_bound(blocks, radii, forcing)
        self.assertLessEqual(lower, upper*(1+1e-12))
        self.assertGreater(lower, 0.)

    def test_attained_bound_dominates_every_sampled_admissible_state(self):
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
        blocks, radii, forcing = self.example()
        shares = D.source_contributions(blocks, radii, forcing)
        self.assertEqual(shares["forcing_template"], float(np.linalg.norm(forcing)))
        self.assertAlmostEqual(sum(shares.values()),
                               D.certified_upper_bound(blocks, radii, forcing), places=12)

    def test_retention_is_positively_homogeneous_in_the_declared_radii(self):
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
        radii = D.declared_radii()
        blocks = {"attitude": np.eye(3)}
        active = {"attitude": radii["attitude"]}
        forcing = np.zeros(3)
        attained = D.attained_lower_bound(blocks, active, forcing)
        self.assertLessEqual(attained/radii["attitude"], 1.+D.MARGIN)
        self.assertGreater(D.MARGIN, 0.)

    def test_restricting_the_initial_set_can_only_lower_the_bound(self):
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
        transitions, responses, covariance = self.ellipsoid_example()
        radii = D.declared_radii()
        base = D.ellipsoid_retention(transitions, responses, radii, covariance)
        scaled = D.ellipsoid_retention(transitions, responses, radii, 4*covariance)
        for name, item in base["groups"].items():
            self.assertAlmostEqual(scaled["groups"][name]["critical_initial_sigma_level"],
                                   item["critical_initial_sigma_level"]/2., places=9)

    def test_a_group_the_forcing_alone_evicts_admits_no_level(self):
        transitions, responses, covariance = self.ellipsoid_example()
        radii = D.declared_radii()
        evicted = [np.zeros(21), np.zeros(21)]
        evicted[1][3:6] = 2*radii["gyro_bias"]
        out = D.ellipsoid_retention(transitions, evicted, radii, covariance)
        self.assertEqual(out["groups"]["gyro_bias"]["critical_initial_sigma_level"], 0.)
        self.assertEqual(out["limiting_group"], "gyro_bias")
        self.assertFalse(out["P4_PASS"])

    def test_a_zero_gain_prefix_never_limits_the_level(self):
        radii = D.declared_radii()
        covariance = np.eye(21)
        transitions = [np.zeros((21, 21)), np.eye(21)]
        responses = [np.zeros(21), np.zeros(21)]
        out = D.ellipsoid_retention(transitions, responses, radii, covariance)
        for name, item in out["groups"].items():
            self.assertAlmostEqual(item["critical_initial_sigma_level"],
                                   radii[name], places=9)


if __name__ == "__main__":
    unittest.main()
