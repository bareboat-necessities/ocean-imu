import unittest
from tools.stability.ou3_theorem.planar_service_enclosure import candidate,oracle_service_floor,scheduler_phase_invariant
class PlanarServiceEnclosureTests(unittest.TestCase):
 def test_candidate_has_feasibility_margin_but_fails_closed(self):
  c=candidate();self.assertFalse(c["verified"]);self.assertFalse(c["theorem_closed"])
  self.assertGreater(c["fresh_tail_observations"]["all_sample_root_tail_service_min"],7)
  self.assertEqual(c["fresh_tail_observations"]["parity_off_fro_max"],0)
  self.assertGreater(c["oracle_feasibility"]["magnetic_information_floor"],1)
  self.assertTrue(c["guard_feasibility"]["fundamental_below_engagement"])
  self.assertFalse(c["guard_feasibility"]["startup_transient_and_harmonics_enclosed"])
  self.assertGreater(len(c["open_dependencies"]),0)
  self.assertTrue(c["scheduler_phase_coordinate_invariant"])
  self.assertTrue(scheduler_phase_invariant(.05,.1))
  self.assertFalse(scheduler_phase_invariant(.1,.1))
 def test_oracle_margin_is_not_promoted(self):
  self.assertGreater(oracle_service_floor(),1)
if __name__=="__main__":unittest.main()
