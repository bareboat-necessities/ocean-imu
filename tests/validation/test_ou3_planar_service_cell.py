import unittest
import numpy as np
from tools.stability.ou3_theorem.planar_service_cell import (
    linked_shift, linked_relative_cell_bound, aw_floor, certificate)
from tools.stability.ou3_theorem.scheduler_s_shift_bound import defect, norm_bound
from tools.stability.ou3_theorem.planar_service_frechet import symmetric_basis, apply_maps
from tools.stability.ou3_theorem.planar_periodic_covariance_tube import certificate as old_tube

class PlanarServiceCellTests(unittest.TestCase):
    def test_exact_certificates_fail_closed_on_admission(self):
        c=certificate()
        self.assertEqual(c['linked_identity_exact_rational_residual'],'0')
        self.assertFalse(c['all_time_magnetic_service_verified'])
        self.assertFalse(c['scheduler_phase_cell_forward_invariant'])
        self.assertFalse(c['aw_floor']['loewner_monotone'])

    def test_old_covariance_order_bound_is_really_false(self):
        P=np.array([[5000.,2.5],[2.5,.005]])
        upper=np.diag([10000.,.01]);Q=np.diag([0.,.01]);H=np.array([[0.,1.]])
        R=np.array([[.01]]);A=np.eye(2)
        D=defect(P,A,Q,H,R)
        self.assertAlmostEqual(D[0,0],500/3,places=8)
        old=(np.linalg.norm(A@upper@H.T,2)**2+np.linalg.norm((upper+Q)@H.T,2)**2)/.01
        self.assertLess(old,1.)
        self.assertGreaterEqual(norm_bound(P,upper,A,Q,H,R),np.linalg.norm(D,2))

    def test_linked_shift_and_rank_with_noncommuting_innovations(self):
        rng=np.random.default_rng(7401)
        for n,m in [(12,2),(9,1),(5,2)]:
            X=rng.normal(size=(n,n));P=X@X.T+np.eye(n)
            A=np.eye(n)+.005*rng.normal(size=(n,n));Q=.003*np.eye(n)
            H=rng.normal(size=(m,n));R=np.eye(m)
            D,info=linked_shift(P,A,Q,H,R)
            np.testing.assert_allclose(D,defect(P,A,Q,H,R),rtol=1e-10,atol=2e-13)
            self.assertLessEqual(np.linalg.matrix_rank(D,tol=1e-10),2*m)
            self.assertAlmostEqual(info['norm_diagnostic'],np.linalg.norm(D,2),places=12)

    def test_linked_cell_bound_retains_one_covariance(self):
        rng=np.random.default_rng(7453);n=5;m=2
        X=rng.normal(size=(n,n));P=X@X.T+np.eye(n);L=np.linalg.cholesky(P)
        A=np.eye(n)+.005*rng.normal(size=(n,n));Q=.002*np.eye(n)
        H=rng.normal(size=(m,n));R=np.eye(m);W=np.linalg.inv(L);radius=.03
        bound=linked_relative_cell_bound(P,L,radius,A,Q,H,R,W)
        for _ in range(20):
            E=rng.normal(size=(n,n));E=(E+E.T)/2;E*=radius/np.linalg.norm(E,2)
            D,_=linked_shift(P+L@E@L.T,A,Q,H,R,W)
            self.assertLessEqual(np.linalg.norm(D,2),bound)
        with self.assertRaises(ValueError):linked_relative_cell_bound(P,L,2,A,Q,H,R)

    def test_aw_floor_is_not_loewner_monotone(self):
        P=np.eye(2);U=P+np.ones((2,2));T=np.array([[3.]])
        gap=aw_floor(U,T,[1])-aw_floor(P,T,[1])
        np.testing.assert_array_equal(gap,[[1.,1.],[1.,0.]])
        self.assertLess(np.linalg.det(gap),0.)

    def test_aw_floor_frobenius_and_target_bounds(self):
        rng=np.random.default_rng(8406)
        for _ in range(30):
            X=rng.normal(size=(6,6));P=X@X.T
            Y=rng.normal(size=(6,6));U=Y@Y.T
            T=np.diag(rng.uniform(.1,4,2));S=np.diag(rng.uniform(.1,4,2))
            fixed=np.linalg.norm(aw_floor(P,T,[3,4])-aw_floor(U,T,[3,4]),'fro')
            self.assertLessEqual(fixed,np.linalg.norm(P-U,'fro')+1e-12)
            varied=np.linalg.norm(aw_floor(P,T,[3,4])-aw_floor(U,S,[3,4]),'fro')
            self.assertLessEqual(varied,np.linalg.norm(P-U,'fro')+np.linalg.norm(T-S,'fro')+1e-12)

    def test_symmetric_basis_is_orthonormal(self):
        B=symmetric_basis(5)
        np.testing.assert_allclose(np.einsum('aij,bij->ab',B,B),np.eye(15),atol=3e-16)

    def test_aw_derivative_drops_only_aw_marginal(self):
        from tools.stability.ou3_theorem.planar_parity import ODD
        n=len(ODD);X=np.ones((1,n,n));idx=ODD.index(16)
        Y=apply_maps(X,[(np.eye(21),True)],ODD)
        expected=X.copy();expected[0,idx,idx]=0
        np.testing.assert_array_equal(Y,expected)
        self.assertEqual(Y[0,idx,0],1.)

    def test_old_tube_does_not_claim_missing_shipping_maps(self):
        c=old_tube()
        self.assertFalse(c['literal_reset_AW_and_joint_coefficient_stream_attached'])
        self.assertFalse(c['AW_floor_Loewner_monotonicity'])
        self.assertIsNone(c['required_retention'])

if __name__=='__main__':unittest.main()
