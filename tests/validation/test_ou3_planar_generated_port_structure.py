import unittest
from tools.stability.ou3_theorem.planar_generated_port_structure import certificate
class GeneratedPortStructureTest(unittest.TestCase):
 def test_default_tuner_has_no_mekf_feedback(self):
  c=certificate()
  self.assertEqual(c["mekf_mean_or_covariance_to_tuner_gain"],0)
  self.assertEqual(c["mekf_mean_or_covariance_to_S_period_gain"],0)
  self.assertEqual(c["mekf_mean_or_covariance_to_AW_target_gain"],0)
  self.assertTrue(c["tuner_to_covariance_is_one_way"])
  self.assertFalse(c["all_time_magnetic_service_verified"])
 def test_locked_live_reference_port_has_its_own_qualified_source_proof(self):
  c=certificate()["locked_live_reference"]
  self.assertEqual(c["locked_Live_reference_reverse_MEKF_port"],0)
  self.assertEqual(c["locked_Live_corrected_mag_reverse_MEKF_port"],0)
  self.assertTrue(c["live_gravity_gate_reads_MEKF"])
  self.assertFalse(c["live_gravity_gate_controls_locked_refinement_or_slew"])
  self.assertTrue(c["refinement_yaw_write_remains_actual_estimator_jump"])
  self.assertFalse(c["source_reference_derivatives_discarded"])
  self.assertFalse(c["uniform_complete_net_work_domination"])
if __name__=="__main__":unittest.main()
