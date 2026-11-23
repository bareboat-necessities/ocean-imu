import unittest
from unittest.mock import patch
from tools.stability.ou3_theorem import planar_generated_port_structure as ports
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
  self.assertEqual(c["real_wave_handover_reference_norm_over_delivered_field"],"101/100")
  self.assertEqual(c["locked_Live_reference_reverse_MEKF_port"],0)
  self.assertEqual(c["locked_Live_corrected_mag_reverse_MEKF_port"],0)
  self.assertEqual(c["live_gravity_gate_MEKF_read_condition"],"with_mag && !mag_ref_set_")
  self.assertFalse(c["locked_Live_gravity_gate_reads_MEKF"])
  self.assertFalse(c["live_gravity_gate_controls_locked_refinement_or_slew"])
  self.assertTrue(c["refinement_yaw_write_remains_actual_estimator_jump"])
  self.assertFalse(c["source_reference_derivatives_discarded"])
  self.assertFalse(c["uniform_complete_net_work_domination"])
 def test_locked_reference_audit_rejects_unconditional_gravity_gate(self):
  shipping=ports.PROXY.read_text()
  with patch.object(ports, "PROXY") as proxy:
   proxy.read_text.return_value=shipping.replace("cfg_.with_mag && !mag_ref_set_", "cfg_.with_mag")
   with self.assertRaisesRegex(ValueError,"gravity-gate"):
    ports.locked_live_reference_audit()
if __name__=="__main__":unittest.main()
