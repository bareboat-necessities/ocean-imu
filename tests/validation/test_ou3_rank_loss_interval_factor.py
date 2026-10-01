from __future__ import annotations
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.rank_loss_interval_factor import (
    FactorState,exact,eye,zeros,verified_inverse,schur_scalar_information,residualized_gram,generalized_ratio_lower,
    source_range_audit,
)

class RankLossIntervalFactorTests(unittest.TestCase):
    def test_literal_four_s_point_box_certifies(self):
        from tools.stability.ou3_theorem.rank_loss_literal_boxes import four_s_gamma_box
        cells=((.10,.10),(.35,.35),(.65,.65),(.95,.95))
        z=four_s_gamma_box(cells,(1.0,1.0))
        self.assertTrue(z["verified"])
        self.assertGreater(z["lower"],0.0)

    def test_linked_soft_return_uses_same_direction(self):
        from tools.stability.ou3_theorem.linked_soft_return import linked_cauchy,rank_one_soft_return
        self.assertTrue(linked_cauchy(2.0,3.0,2.0)["cauchy_holds"])
        z=rank_one_soft_return(1.0,4.0,.5,.2,3.0)
        self.assertGreater(z["effective_information"],0.0)
        self.assertLess(z["posterior_rank_one_variance"],z["propagated_prior_variance"])

    def test_adjacent_kernel_line_recurrence(self):
        from tools.stability.ou3_theorem.adjacent_kernel_lines import adjacent_line_update,kernel_carry_derivative
        a=4.0;c=3.0;j=.2
        ortho=adjacent_line_update(a,0.0,j,c)
        mid=adjacent_line_update(a,.6,j,c)
        same=adjacent_line_update(a,1.0,j,c)
        self.assertEqual(ortho["next_kernel_variance"],0.0)
        self.assertLess(mid["next_kernel_variance"],same["next_kernel_variance"])
        self.assertTrue(same["kernel_invariance"])
        self.assertGreater(kernel_carry_derivative(a,.4,j,c),0.0)

    def test_scalar_kernel_regularization_removes_J_zero_blowup(self):
        from tools.stability.ou3_theorem.kernel_regularized_scalar import regularized_scalar,augmented_rank_one_leverage
        z=regularized_scalar(0.0,4.0,3.0)
        self.assertAlmostEqual(z["terminal_excess"],12.0)
        self.assertAlmostEqual(z["upper_cB"],12.0)
        a=augmented_rank_one_leverage(4.0,3.0)
        self.assertLess(a["leverage"],1.0)
        self.assertEqual(augmented_rank_one_leverage(4.0,None)["leverage_upper"],1.0)

    def test_terminal_compatibility_energy_composition(self):
        from tools.stability.ou3_theorem.terminal_compatibility_energy import scalar_terminal_bound,source_range_status
        z=scalar_terminal_bound(2.0,3.0,4.0,5.0)
        self.assertAlmostEqual(z["B_term"],2.6)
        self.assertAlmostEqual(z["C_c_upper"],10.4)
        self.assertFalse(source_range_status(False)["verified"])
        self.assertTrue(source_range_status(True,7.0)["verified"])

    def test_homogeneous_compatibility_scalar_ratio_scale_invariant(self):
        from tools.stability.ou3_theorem.compatibility_chart import reduced_pair,scalar_ratio,kernel_prior_floor_for_unit_direction
        R=exact([[2.0],[0.0]])
        Qnm=zeros(2,2);N=exact([[3.0,0],[0,0]])
        # Projector onto span(e1).
        Pk=exact([[1.0,0],[0,0]])
        q,n=reduced_pair(R,Qnm,N,Pk,4.0)
        z=scalar_ratio(q.mid[0][0],n.mid[0][0])
        self.assertTrue(z["verified"]);self.assertAlmostEqual(z["upper"],12.0)
        self.assertAlmostEqual(kernel_prior_floor_for_unit_direction(4.0),.25)

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


class LinkedSoftFullCovarianceTests(unittest.TestCase):
    """Exact algebraic witnesses; none is a shipping-reachability certificate."""

    def setUp(self):
        from fractions import Fraction
        from tools.stability.ou3_theorem.lin_path_certificate import inverse
        from tools.stability.ou3_theorem.matrix_certificates import add, identity, matmul, transpose
        self.F = Fraction
        self.inverse, self.add = inverse, add
        self.eye, self.mm, self.tr = identity, matmul, transpose

    def matrix(self, rows):
        return [[self.F(x) for x in row] for row in rows]

    def scale(self, rows, scalar):
        return [[scalar*x for x in row] for row in rows]

    def quad(self, matrix, vector):
        v = [[self.F(x)] for x in vector]
        return self.mm(self.tr(v), self.mm(matrix, v))[0][0]

    def test_full_baseline_rank_one_identity_with_singular_and_full_information(self):
        from tools.stability.ou3_theorem.linked_soft_return import conditional_rank_one
        f = self.F
        b = self.matrix([[2, '1/3'], ['1/3', 3]])
        y = [f(3, 5), f(4, 5)]
        yy = self.mm([[x] for x in y], [y])
        for j in (self.matrix([[0, 0], [0, 2]]), self.matrix([[1, '1/4'], ['1/4', 2]])):
            for a in (f(0), f(7, 3)):
                z = conditional_rank_one(b, j, y, a)
                prior = self.add(b, self.scale(yy, a))
                direct = self.inverse(self.add(self.inverse(prior), j))
                self.assertEqual(z['posterior'], direct)
                self.assertGreaterEqual(z['effective_information'], 0)

    def test_sharp_generalized_eigenvalue_is_attained(self):
        from tools.stability.ou3_theorem.linked_soft_return import linked_directional_max
        f = self.F
        b = self.matrix([[2, '1/3'], ['1/3', 3]])
        j = self.matrix([[1, '1/4'], ['1/4', 2]])
        n, p = [f(5, 13), f(12, 13)], f(7, 3)
        for e in (None, self.matrix([[2, 1], [0, 3]]), self.matrix([[3], [4]])):
            z = linked_directional_max(b, j, n, p, e)
            extremizer = z['extremizing_prior']
            actual = self.inverse(self.add(self.inverse(extremizer), j))
            self.assertEqual(self.quad(actual, n), z['max_variance'])
            excess = self.add(extremizer, b, -1)
            self.assertEqual(excess[0][0]+excess[1][1], p)
            self.assertEqual(excess[0][0]*excess[1][1]-excess[0][1]**2, 0)

    def test_full_space_duality_equals_isotropic_envelope_directionally(self):
        from tools.stability.ou3_theorem.linked_soft_return import linked_directional_max
        b = self.matrix([[1, '-9/10'], ['-9/10', 1]])
        j = self.matrix([[0, 0], [0, 1]])
        p = self.F(3, 2)
        z = linked_directional_max(b, j, [1, 0], p)
        direct = self.inverse(self.add(self.inverse(self.add(b, self.scale(self.eye(2), p))), j))
        self.assertEqual(z['max_variance'], direct[0][0])

    def test_nonorthonormal_chart_scaling_preserves_sharp_bound(self):
        from tools.stability.ou3_theorem.linked_soft_return import linked_directional_max
        b, j = self.eye(2), self.matrix([[0, 0], [0, 2]])
        values = [linked_directional_max(b, j, [1, 1], 3, self.matrix([[3*k], [4*k]]))
                  for k in (self.F(1), self.F(1, 1000), self.F(1000))]
        self.assertEqual(values[0]['max_variance'], values[1]['max_variance'])
        self.assertEqual(values[0]['max_variance'], values[2]['max_variance'])

    def test_correlated_baseline_changes_worst_angle(self):
        from tools.stability.ou3_theorem.linked_soft_return import conditional_rank_one, linked_directional_max
        f = self.F
        b, j = self.matrix([[1, '-9/10'], ['-9/10', 1]]), self.matrix([[0, 0], [0, 1]])
        aligned = conditional_rank_one(b, j, [1, 0], 1)
        rotated = conditional_rank_one(b, j, [f(12, 13), f(5, 13)], 1)
        extra = rotated['posterior'][0][0]-rotated['baseline'][0][0]
        self.assertEqual(extra, f(1083, 968))
        self.assertGreater(rotated['posterior'][0][0], aligned['posterior'][0][0])
        bound = linked_directional_max(b, j, [1, 0], 1)
        self.assertEqual(bound['max_excess'], f(227, 200))
        self.assertEqual(bound['max_variance'], f(173, 100))
        self.assertLessEqual(rotated['posterior'][0][0], bound['max_variance'])

    def test_compressed_covariance_is_not_compressed_inverse(self):
        pi = self.matrix([[1, '4/5'], ['4/5', 1]])
        linked_product = pi[0][0]*self.inverse(pi)[0][0]
        self.assertEqual(linked_product, self.F(25, 9))
        self.assertGreater(linked_product, 1)  # The one-dimensional compression has condition number 1.

    def test_product_can_fail_while_exact_scalar_return_is_invariant(self):
        from tools.stability.ou3_theorem.linked_soft_return import soft_scalar_return
        # Pi=diag(100,1), Phi=I/2, J=diag(0,1), n=u=e1.
        z = soft_scalar_return(100, self.F(1, 4), 0, 200)
        self.assertEqual(z['variance'], 150)
        self.assertEqual(z['ratio'], self.F(3, 4))
        self.assertGreater(self.F(100, 4), 1)  # d_soft*H=25 on this very same word.

    def test_synthetic_next_precision_does_not_prove_invariance(self):
        from tools.stability.ou3_theorem.linked_soft_return import (
            conditional_rank_one, rank_one_soft_return, soft_scalar_return)
        b, j = self.eye(2), self.matrix([[0, 0], [0, 1]])
        for c in (1, 2, 100):
            real = conditional_rank_one(b, j, [1, 0], c)
            self.assertEqual(real['effective_information'], 0)
            self.assertEqual(real['posterior'][0][0], c+1)
            exact = soft_scalar_return(1, 1, 0, c)
            self.assertEqual(exact['invariance_margin'], -1)
            synthetic = rank_one_soft_return(c, 1, 1, 1, c)
            self.assertEqual(synthetic['next_kernel_variance'], c/2)
            self.assertFalse(synthetic['proves_shipping_invariance'])

    def test_exact_spectral_split_uses_full_quotient_and_baseline(self):
        from tools.stability.ou3_theorem.linked_soft_return import soft_scalar_return
        f = self.F
        u, q, n = [f(3, 5), f(4, 5)], [f(-4, 5), f(3, 5)], [f(5, 13), f(12, 13)]
        uu, qq = self.mm([[x] for x in u], [u]), self.mm([[x] for x in q], [q])
        pi, phi, c = self.matrix([[2, '1/3'], ['1/3', 3]]), self.matrix([[2, 1], [-1, 3]]), f(7, 5)
        w = [row[0] for row in self.mm(self.tr(phi), [[x] for x in n])]
        d = self.quad(pi, n)+sum(x*y for x, y in zip(w, q))**2/2
        ell_sq = sum(x*y for x, y in zip(w, u))**2
        for soft_info in (f(0), f(3, 11)):
            j = self.add(self.scale(uu, soft_info), self.scale(qq, 2))
            g = self.add(j, self.scale(uu, 1/c))
            direct = self.add(pi, self.mm(phi, self.mm(self.inverse(g), self.tr(phi))))
            self.assertEqual(soft_scalar_return(d, ell_sq, soft_info, c)['variance'], self.quad(direct, n))

    def test_physical_prior_scale_and_polynomial_identity(self):
        from tools.stability.ou3_theorem.linked_soft_return import soft_scalar_return
        f = self.F
        for j in (f(0), f(1, 7)):
            z = soft_scalar_return(f(2, 3), f(5, 4), j, 3, 4)
            self.assertEqual(z['fixed_point_polynomial'], z['invariance_margin']*(4+3*j))
        self.assertGreater(soft_scalar_return(1, 4, 1, 5)['invariance_margin'], 0)
        self.assertLess(soft_scalar_return(1, 4, 1, 4)['invariance_margin'], 0)

    def test_invalid_domains_fail_closed(self):
        from tools.stability.ou3_theorem.linked_soft_return import conditional_rank_one, linked_directional_max, soft_scalar_return
        with self.assertRaises(ValueError):
            conditional_rank_one(self.eye(2), [[0, 0], [0, -1]], [1, 0], 1)
        with self.assertRaises(ValueError):
            linked_directional_max(self.eye(2), self.eye(2), [1, 0], 1, [[1, 1], [0, 0]])
        for c in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                soft_scalar_return(1, 1, 0, c)

    def test_status_never_promotes_algebra_to_O2(self):
        from tools.stability.ou3_theorem.linked_soft_return import linked_product_status
        z = linked_product_status()
        self.assertTrue(z['full_covariance_baseline_required'])
        self.assertFalse(z['next_ceiling_is_an_applied_measurement'])
        self.assertFalse(z['source_uniform_O2_closed'])


if __name__=="__main__":unittest.main()
