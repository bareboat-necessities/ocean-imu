"""Universal polynomial identities and finite runtime chart regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from test_finite_measurement_graph import Poly
from test_finite_core import root, physical_successor
from test_finite_fresh_joint24_entry import entry, reference
from tools.stability.ou3_alt_contraction import finite_attitude_atlas as A
from tools.stability.ou3_alt_contraction import finite_core as C
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE


def represented_state(chart, mode='H'):
    s = root(mode)
    q = A.homogeneous((F(1,10), F(-1,20), F(1,30)), chart)
    z = (F(1,10), F(-1,20), F(1,30))+s.z[3:]
    return replace(s, z=z, attitude_chart=chart,
                   reference=replace(s.reference, q_world_to_body=q))


class Tests(unittest.TestCase):
    def test_rotation_orthogonality_and_nonzero_product_universal_polynomials(self):
        with patch.object(Poly, 'N', 8):
            q = [Poly.var(i) for i in range(4)]
            r = [Poly.var(i) for i in range(4,8)]
            n2 = M.dot(q,q)
            N = A.rotation_numerator(q)
            self.assertEqual(M.mm(M.transpose(N),N), M.scaled(M.eye(3),n2*n2))
            qr = P.quat_mul(q,r)
            self.assertEqual(M.dot(qr,qr), n2*M.dot(r,r))
            self.assertEqual(A.rotation_numerator(qr), M.mm(N,A.rotation_numerator(r)))

    def test_all_chart_transport_coefficients_all_indeterminates(self):
        with patch.object(Poly, 'N', 11):
            c = [Poly.var(i) for i in range(3)]
            left = [Poly.var(i) for i in range(3,7)]
            right = [Poly.var(i) for i in range(7,11)]
            for chart in range(4):
                q = [Poly(2) if i == chart else Poly(0) for i in range(4)]
                for j, i in enumerate(A.slots(chart)): q[i] = c[j]
                qp = P.quat_mul(P.quat_mul(left,q),right)
                columns = [P.quat_mul(P.quat_mul(left,[F(i == j) for i in range(4)]),right)
                           for j in range(4)]
                matrix = M.transpose(columns)
                for target in range(4):
                    for i in A.slots(target):
                        numerator = 4*matrix[i][chart]+sum(
                            2*matrix[i][j]*x for j,x in zip(A.slots(chart),c))
                        self.assertEqual(numerator,2*qp[i])

    def test_cover_exact_poles_nearby_family_and_projective_scaling(self):
        family = [(0,0,0,1), (0,1,0,0), (0,0,1,0), (1,1,-1,-1)]
        family += [(F(2*n,n*n+1),0,0,F(n*n-1,n*n+1)) for n in (1,2,10,10**30)]
        for q in family:
            p = A.encode(q)
            self.assertLessEqual(max(abs(x) for x in p.coordinates),2)
            self.assertGreaterEqual(4*q[p.chart]**2, M.dot(q,q))
            self.assertEqual(A.rotation(p),C.rotation(q))
            for scale in (F(-7,3),F(1,10**60)):
                self.assertEqual(A.encode([scale*x for x in q]),p)
        self.assertEqual(A.encode((1,1,1,1)).chart,0)
        with self.assertRaises(ValueError): A.encode((0,0,0,0))
        with self.assertRaises(ValueError): A.encode((0,0,0,1),chart=0)

    def test_every_overlap_and_affine_transport_preserve_rotation(self):
        q=(1,2,-3,4)
        for source in range(4):
            p=A.encode(q,chart=source)
            for target in range(4):
                tr=A.transport(p,chart_after=target)
                self.assertEqual(tr.after,A.encode(q,chart=target))
                self.assertEqual(A.rotation(tr.before),A.rotation(tr.after))
                self.assertEqual(tuple(M.mv(tr.affine,(*p.coordinates,1))),tr.after.coordinates)
        p=A.encode((0,0,0,1))
        tr=A.transport(p,right=(0,0,0,-1))
        self.assertEqual(tr.after,A.Point(0,(0,0,0)))

    def test_fresh_entry_covers_admitted_south_without_changing_any_other_slot(self):
        e=entry(); r=replace(reference(),q_world_to_body=(0,0,0,1))
        south=FRESH.build(e,r,scope=SCOPE.certified_scope())
        north=FRESH.build(e,reference(),scope=SCOPE.certified_scope())
        self.assertEqual(south.attitude_chart,3)
        self.assertEqual(south.z[:3],(0,0,0))
        self.assertEqual(south.z[3:],north.z[3:])
        self.assertEqual(south.covariance,e.P)
        self.assertEqual(south.q_hat,e.q_hat)
        self.assertIs(south.reference,r)
        # Equal numerical coordinate zero in these two charts is NOT equal error.
        self.assertNotEqual(A.rotation(A.Point(3,south.z[:3])),M.eye(3))

    def test_prediction_all_charts_and_bias_modes_keep_same_physical_recurrence(self):
        for chart in range(4):
            for mode in ('H','A'):
                s=represented_state(chart,mode)
                seg,kw=physical_successor(s)
                out=C.prediction(s,seg,**kw)
                expected=A.encode(P.quat_mul(seg.after.q_world_to_body,P.quat_conj(out.q_hat)))
                self.assertEqual((out.attitude_chart,out.z[:3]),(expected.chart,expected.coordinates))
                self.assertEqual(out.z[21:],seg.after.beta)
                self.assertIs(out.reference,seg.after)
                self.assertEqual(out.covariance,s.covariance)

    def test_all_sensor_kinds_charts_and_modes_retain_actual_solve_and_nominal_update(self):
        for chart in range(4):
            for mode in ('H','A'):
                s=represented_state(chart,mode)
                for kind,packet in (
                    ('S_zero',{}),
                    ('accelerometer',{'observed':(F(11,100),0,F(-196133,20000))}),
                    ('magnetometer',{'observed':(22,F(1,10000),43),'magnetic_reference':(22,0,43)})):
                    result=C.measurement(s,kind,R=M.eye(3),**packet)
                    out=result.state
                    expected=A.encode(P.quat_mul(s.reference.q_world_to_body,P.quat_conj(out.q_hat)))
                    self.assertEqual((out.attitude_chart,out.z[:3]),(expected.chart,expected.coordinates))
                    self.assertEqual(out.z[3:21],tuple(s.z[i]-result.correction[i] for i in range(3,21)))
                    self.assertEqual(out.z[21:],s.z[21:])
                    self.assertEqual(M.mm(result.gain,result.innovation),list(map(list,result.numerator)))

    def test_H18_A21_hold_covariance_edges_preserve_nonzero_chart(self):
        s=represented_state(3)
        cfg=GATE.Config()
        control=GATE.State(250,F(1),False,True)
        release=GATE.set_hold(control,s,cfg,hold=False,live=True)
        held=GATE.set_hold(release.state,release.filter_state,cfg,hold=True,live=True)
        for result in (release,held):
            self.assertEqual(result.filter_state.attitude_chart,3)
            self.assertEqual(result.filter_state.z,s.z)
            self.assertEqual(result.filter_state.q_hat,s.q_hat)

    def test_chart_tampering_is_rejected_by_existing_core(self):
        s=represented_state(3)
        with self.assertRaises(ValueError): replace(s,attitude_chart=0)
        with self.assertRaises(ValueError): replace(s,attitude_chart=True)
        with self.assertRaises(ValueError): replace(s,attitude_chart=4)


if __name__=='__main__': unittest.main()
