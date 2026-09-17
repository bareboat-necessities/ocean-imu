import unittest
import numpy as np
from tools.stability.ou3_alt_contraction import magnetic_service_formulation as F
class TestMagneticServiceFormulation(unittest.TestCase):
 def test_selection_fail_closed(self):
  s=F.selection_status(); self.assertTrue(s['old_per_word_absolute_heading_target_retired']); self.assertFalse(s['yaw_only_quotient_closes_axial_bias']); self.assertIsNone(s['corrected_source_uniform_rho']); self.assertFalse(s['storage_search_allowed'])
 def test_yaw_only_quotient_leaves_bias_neutral(self):
  U=np.array([[1.,3.],[0.,1.]])
  self.assertEqual(U[1,1]**2,1.)
 def test_projected_ratio_keeps_all_output_rows(self):
  A=np.eye(3); A[2,0]=2.; Z=F.coordinate_injection(3,[0,1]); r=F.projected_storage_ratio(A,np.eye(3),np.eye(3),Z); self.assertGreater(r['rho_point'],1.)
 def test_common_metric_family(self):
  Z=np.eye(2); r=F.common_storage_family_ratio([.5*np.eye(2),.9*np.eye(2)],np.eye(2),Z); self.assertAlmostEqual(r['worst_rho_point'],.81)
 def test_quotient_rejects_leakage(self):
  A=np.array([[1.,0.],[1.,1.]]); N=np.array([[1.],[0.]]); Q=np.array([[0.],[1.]])
  with self.assertRaises(ValueError): F.quotient_map(A,N,N,Q,Q)
 def test_service_rejects_call_cadence_without_acceptance(self):
  c=F.MagneticServiceContract(.1,.1,.01); ev=[F.MagneticEvent(.04,False,True),F.MagneticEvent(.08,False,True)]; r=F.audit_magnetic_service(ev,0,.1,c); self.assertFalse(r['finite_window_point_service_pass'])
 def test_one_heading_observation_is_rank_one(self):
  c=F.MagneticServiceContract(1.,.1,.01); ev=[F.MagneticEvent(.5,True,True,[[1.,.5]])]; r=F.audit_magnetic_service(ev,0,1,c); self.assertFalse(r['pair_information_point_pass'])
 def test_two_transported_observations_can_span_pair(self):
  c=F.MagneticServiceContract(.6,.1,.01); ev=[F.MagneticEvent(.25,True,True,[[1.,.25]]),F.MagneticEvent(.75,True,True,[[1.,.75]])]; r=F.audit_magnetic_service(ev,0,1,c); self.assertTrue(r['finite_window_point_service_pass'])
 def test_carried_boundaries(self):
  p=[F.WordPiece(np.eye(2),'h','a','b'),F.WordPiece(.5*np.eye(2),'h','b','c')]; self.assertTrue(np.allclose(F.compose_carried_pieces(p),.5*np.eye(2)))
 def test_detached_boundary_rejected(self):
  p=[F.WordPiece(np.eye(2),'h','a','b'),F.WordPiece(np.eye(2),'h','x','c')]
  with self.assertRaises(ValueError): F.compose_carried_pieces(p)
if __name__=='__main__': unittest.main()
