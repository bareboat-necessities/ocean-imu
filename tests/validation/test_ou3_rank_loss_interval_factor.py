from __future__ import annotations
import math,sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.rank_loss_interval_factor import (
    FactorState,exact,eye,zeros,verified_inverse,source_covariance,
    schur_scalar_information,residualized_gram,generalized_ratio_lower,
    source_range_audit,
)

class RankLossIntervalFactorTests(unittest.TestCase):
    def test_literal_four_s_point_box_certifies(self):
        from tools.stability.ou3_theorem.rank_loss_literal_boxes import four_s_gamma_box
        cells=((.10,.10),(.35,.35),(.65,.65),(.95,.95))
        z=four_s_gamma_box(cells,(1.0,1.0))
        self.assertTrue(z["verified"])
        self.assertGreater(z["lower"],0.0)

    def test_complete_magnetic_cancellation_reduction(self):
        from tools.stability.ou3_theorem.complete_magnetic_gram import scalar_magnetic_cancellation,restricted_complete_ratio
        z=scalar_magnetic_cancellation(2.0,3.0)
        self.assertTrue(z["cancellable"]);self.assertLess(z["selected_direction_action"],1e-28)
        R=eye(2);q=exact([[2,0],[0,3]]);n=exact([[1,0],[0,1]])
        r=restricted_complete_ratio(R,q,n)
        self.assertTrue(r["verified"]);self.assertGreater(r["lower"],1.9)

    def test_relative_process_modulus_has_zero_uniform_limit(self):
        from tools.stability.ou3_theorem.magnetic_service_schur import relative_process_limit,relative_process_route_status
        a=relative_process_limit(1.0,1.0);b=relative_process_limit(1e12,1.0)
        self.assertAlmostEqual(a["joint_gamma"],.5)
        self.assertLess(b["joint_gamma"],1e-11)
        self.assertEqual(relative_process_route_status()["source_uniform_q_rel_lower"],0.0)

    def test_aggregate_magnetic_service_counterexample(self):
        from tools.stability.ou3_theorem.magnetic_service_schur import service_only_counterexample,joint_process_modulus
        z=service_only_counterexample()
        self.assertEqual(z["service_mu"],1.0)
        self.assertEqual(z["canonical_correlation"],1.0)
        self.assertEqual(z["gamma_residualized"],0.0)
        q=joint_process_modulus(1.0)
        self.assertTrue(q["verified"]);self.assertAlmostEqual(q["gamma_joint_lower"],.5)

    def test_theorem_magnetic_seed_audit_fails_on_real_blockers(self):
        from tools.stability.ou3_theorem.magnetic_seed_cover import magnetic_seed_audit,theorem_seed_cover
        z=magnetic_seed_audit()
        self.assertFalse(z["seed_cover_constructible"])
        self.assertIn("full_covariance_P",z["blockers"])
        self.assertIn("event_schedule",z["blockers"])
        seeds,a=theorem_seed_cover()
        self.assertEqual(seeds,())
        self.assertEqual(a["qualification"],"OU3_MAGNETIC_SEED_AUDIT_V1")

    def test_exhaustive_magnetic_strata_driver_positive_exact_leaf(self):
        from tools.stability.ou3_theorem.magnetic_strata_certificate import ClosedMagneticStratum,exhaustive_certificate
        from tools.stability.ou3_theorem.magnetic_literal_box_export import CorrectionBox
        # Three root columns: protected 0,1 and orthogonal nuisance 2.
        H=exact([[1,0,0],[0,1,0],[0,0,1]])
        st=ClosedMagneticStratum("s",3,(CorrectionBox(eye(3),"mag",True,H,eye(3),.5,.5),),
                                 (0,1),(2,),"unit-test-cover")
        z=exhaustive_certificate([st],max_depth=2)
        self.assertTrue(z["verified"]);self.assertGreater(z["gamma_M_lower"],.999999999)

    def test_exhaustive_driver_requires_coverage_tag(self):
        from tools.stability.ou3_theorem.magnetic_strata_certificate import ClosedMagneticStratum,exhaustive_certificate
        from tools.stability.ou3_theorem.magnetic_literal_box_export import CorrectionBox
        H=exact([[1,0],[0,1]])
        st=ClosedMagneticStratum("s",2,(CorrectionBox(eye(2),"mag",True,H,eye(2),.5,.5),),
                                 (0,1),(),"")
        z=exhaustive_certificate([st])
        self.assertFalse(z["verified"]);self.assertEqual(z["gamma_M_lower"],0.0)

    def test_shipping_magnetic_operation_box(self):
        from tools.stability.ou3_theorem.shipping_operation_interval_boxes import magnetic_update_box,reset_box
        P=eye(21); v=exact([[1],[2],[3]]); R=eye(3)
        box,z=magnetic_update_box(P,v,R,applied=True)
        self.assertTrue(z["verified"])
        self.assertEqual(box.H.shape,(3,21));self.assertEqual(box.S_actual.shape,(3,3))
        self.assertEqual(box.A.shape,(21,21))
        d=exact([[.1],[0],[0]])
        self.assertEqual(reset_box(d).G.shape,(21,21))

    def test_shipping_magnetic_update_fails_closed_for_singular_S(self):
        from tools.stability.ou3_theorem.shipping_operation_interval_boxes import magnetic_update_box
        P=zeros(21,21);v=exact([[1],[0],[0]]);R=zeros(3,3)
        box,z=magnetic_update_box(P,v,R,applied=True)
        self.assertFalse(z["verified"])

    def test_literal_magnetic_export_pre_correction_phi(self):
        from tools.stability.ou3_theorem.magnetic_literal_box_export import (
            PredictionBox,CorrectionBox,ResetBox,SyncBox,export_magnetic_event_boxes,event_tuples)
        F=exact([[1,1],[0,1]]); A=exact([[.5,0],[0,1]])
        H=exact([[1,0],[0,1]]); S=eye(2)
        ev=export_magnetic_event_boxes(2,[PredictionBox(F),
            CorrectionBox(A,"mag",True,H,S,.4,.4),ResetBox(eye(2)),SyncBox()])
        self.assertEqual(len(ev),1)
        self.assertEqual(ev[0].Phi_from_window_root.mid,F.mid)
        self.assertEqual(event_tuples(ev)[0][0].mid,H.mid)

    def test_literal_magnetic_rejected_event_not_exported(self):
        from tools.stability.ou3_theorem.magnetic_literal_box_export import (
            CorrectionBox,export_magnetic_event_boxes)
        ev=export_magnetic_event_boxes(2,[CorrectionBox(eye(2),"mag",False,None,None,.2,.2)])
        self.assertEqual(ev,())

    def test_magnetic_nuisance_schur_exact(self):
        from tools.stability.ou3_theorem.magnetic_nuisance_interval import service_schur_information
        # Two protected columns orthogonal to one nuisance column.
        H=exact([[1,0,0],[0,1,0],[0,0,1]])
        z=service_schur_information([(H,eye(3),eye(3))],[0,1],[2])
        self.assertTrue(z["verified"])
        self.assertGreater(z["gamma_M_lower"],.999999999)

    def test_magnetic_service_floor_alone_fails_closed(self):
        from tools.stability.ou3_theorem.magnetic_nuisance_interval import service_floor_only_contract
        z=service_floor_only_contract()
        self.assertFalse(z["verified"])
        self.assertEqual(z["gamma_M_lower"],0.0)

    def test_literal_source_certificate_fails_closed_on_magnetic(self):
        from tools.stability.ou3_theorem.rank_loss_literal_boxes import rank_loss_factor_certificate
        z=rank_loss_factor_certificate(max_depth=2)
        self.assertFalse(z["source_uniform_verified"])
        self.assertEqual(z["gamma_M_lower"],0.0)
        self.assertEqual(z["beta_lower"],0.0)

    def test_verified_generic_inverse(self):
        a=exact([[2.0,.1],[.1,1.0]])
        ai,c=verified_inverse(a)
        self.assertTrue(c["verified"])
        self.assertLess(c["residual_ratio_upper"],1e-12)
        self.assertAlmostEqual(ai.mid[0][0],1/1.99,places=12)

    def test_factor_propagation_preserves_shared_source(self):
        st=FactorState(eye(2),zeros(2,0))
        st=st.predict(exact([[1,1],[0,1]]),exact([[1],[0]]))
        o,c=st.observe(exact([[1,0]]),exact([[2]]))
        self.assertEqual(o.shape,(1,2));self.assertEqual(c.shape,(1,2))
        self.assertAlmostEqual(c.mid[0][0],1.0)
        self.assertAlmostEqual(c.mid[0][1],2.0)

    def test_four_row_schur_scalar_exact_case(self):
        # va is orthogonal to span(V0), R=I => gamma=||va||^2=1.
        v0=exact([[1,0,0],[0,1,0],[0,0,1],[0,0,0]])
        va=exact([[0],[0],[0],[1]])
        z=schur_scalar_information(v0,va,eye(4))
        self.assertTrue(z["verified"]);self.assertGreater(z["lower"],.999999999)

    def test_residualized_magnetic_gram(self):
        obs=exact([[1,0],[0,1],[0,0]])
        nuis=exact([[0],[0],[1]])
        g,z=residualized_gram(obs,nuis,eye(3))
        self.assertTrue(z["verified"]);self.assertGreater(z["lower"],.999999999)
        self.assertEqual(g.shape,(2,2))

    def test_generalized_ratio_lower(self):
        z=generalized_ratio_lower(exact([[2,0],[0,3]]),exact([[4,0],[0,5]]))
        self.assertTrue(z["verified"]);self.assertGreaterEqual(z["lower"],.4-1e-15)

    def test_source_range_audit_fails_closed(self):
        z=source_range_audit()
        self.assertFalse(z["source_uniform_verified"])
        self.assertEqual(z["gamma_S_lower"],0.0)
        self.assertEqual(z["gamma_M_lower"],0.0)
        self.assertEqual(z["beta_lower"],0.0)

if __name__=="__main__":unittest.main()
