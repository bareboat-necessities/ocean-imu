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
            self.assertEqual(scaled[name]["prefix_index"], item["prefix_index"])

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


if __name__ == "__main__":
    unittest.main()
