"""Sampling sharpness, fail-closed admission, and full-vector geometry."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.marine_motion import (
    MarineLimits, continuation_admitted, quiet_water_continuation_certificate,
)
from tools.stability.ou3_theorem.sampling_fidelity import (
    certificate, trapezoidal_weights, vector_information_floor,
)
from tools.stability.ou3_theorem.matrix_certificates import add, is_psd, matmul, transpose
from tools.stability.ou3_theorem.lin_path_certificate import inverse


class SamplingFidelityTests(unittest.TestCase):
    def test_admission_requires_jerk_and_absolute_continuity(self):
        limits=MarineLimits(8.1,5.5,8.8,.611,1100)
        c=quiet_water_continuation_certificate()
        self.assertTrue(continuation_admitted(c,limits))
        for bad in (replace(c,jerk_norm_upper_mps3=None),
                    replace(c,jerk_norm_upper_mps3=101),
                    replace(c,acceleration_locally_absolutely_continuous_certified=False)):
            self.assertFalse(continuation_admitted(bad,limits))
        with self.assertRaises(ValueError):
            replace(c,jerk_norm_upper_mps3=float('nan'))

    def test_nonuniform_weights_and_sharp_tent_defect(self):
        steps=[F(1,250),F(3,500),F(1,200)]
        w=trapezoidal_weights(steps)
        self.assertEqual(sum(w),1)
        self.assertTrue(all(x>0 for x in w))
        # Concatenated negative triangular tents have zero sampled endpoints,
        # slope magnitude J, and exact integral -sum(J*h^2/4). This saturates
        # the bound without discontinuous acceleration or uniform timing.
        j=F(100); total=sum(steps)
        integral=-sum(j*h*h/4 for h in steps)
        error=-integral/total
        self.assertEqual(error,sum(j*h*h/(4*total) for h in steps))
        self.assertLessEqual(error,j*max(steps)/4)
        # The bound cannot be inferred from endpoint changes: these are zero.
        self.assertGreater(error,0)

    def test_full_joint_gram_floor_retains_cross_terms(self):
        # Independent exact PSD check on perturbed gravity/field pairs. The
        # magnetic field is allowed to be predominantly vertical; each single
        # vector factor has rank two, while the joint matrix is positive.
        eps=F(1,12); residual=F(7,75)
        floor=vector_information_floor(eps,residual)
        self.assertGreater(floor,0)
        for sx in (-1,1):
            for sy in (-1,1):
                u=[sx*eps,F(0),F(-1)]
                b=[F(1,5),sy*residual,F(9,10)]
                matrix=[[sum(sum(x*x for x in v)*(i==j)-v[i]*v[j]
                             for v in (u,b))-(floor if i==j else 0)
                         for j in range(3)] for i in range(3)]
                self.assertTrue(is_psd(matrix))
        self.assertEqual(vector_information_floor(F(1,2),residual),0)

    def test_exact_thresholds_and_scope(self):
        c=certificate(); alias=c['stationary_sample_alias']; info=c['joint_vector_information']
        self.assertLess(F(alias['hidden_specific_force_upper_mps2']),F(alias['six_degree_chord_lower_mps2']))
        self.assertEqual(info['magnetometer_residual_total_uT'],'7')
        self.assertGreater(F(info['normalized_joint_gram_eigenvalue_lower']),F(6,100000))
        self.assertGreater(F(c['old_witness_jerk_squared_lower']),F(c['jerk_limit_mps3'])**2)
        self.assertFalse(info['full_21_corrected_loss_certified'])
        self.assertFalse(c['shipping_capture_certified'])

    def test_full_metric_supply_with_carried_cross_covariance(self):
        p=[[F(2),F(3,4)],[F(3,4),F(1)]]
        f=[[F(1),F(1,3)],[F(0),F(4,5)]]
        # One process direction, plus positive sync in the other direction.
        q=F(1,7)
        pred=add(matmul(matmul(f,p),transpose(f)),[[F(1,11),F(0)],[F(0),q]])
        transport=[f[0]+[F(0)],f[1]+[F(1)]]
        out=matmul(matmul(transpose(transport),inverse(pred)),transport)
        ip=inverse(p)
        incoming=[ip[0]+[F(0)],ip[1]+[F(0)],[F(0),F(0),1/q]]
        self.assertTrue(is_psd(add(incoming,out,F(-1))))
        # Finite rotation identity: Q has cos(phi)=3/5 about the first axis.
        rotation=[[F(1),F(0),F(0)],[F(0),F(3,5),-F(4,5)],
                  [F(0),F(4,5),F(3,5)]]
        for v in ([F(2),F(3),F(5)],[F(1,5),F(-1,7),F(9,10)]):
            qv=matmul(rotation,[[x] for x in v])
            displacement=sum((x[0]-y)**2 for x,y in zip(qv,v))
            self.assertEqual(displacement,2*(1-F(3,5))*(v[1]**2+v[2]**2))


if __name__=='__main__':
    unittest.main()
