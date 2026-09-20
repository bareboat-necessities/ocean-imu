"""Exact checks for the continuous-history and all-time service counterexample."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.sampled_capture_obstruction import (
    correct, predict, service_certificate, transition, witness_certificate,
)
from tools.stability.ou3_theorem.matrix_certificates import (
    add, identity, is_psd, matmul, transpose,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse


class SampledCaptureTests(unittest.TestCase):
    def test_exact_distinct_attitudes_produce_same_samples(self):
        c=witness_certificate()
        r=[[F(v) for v in row] for row in c['body_to_world_rotation_plus']]
        self.assertEqual(matmul(transpose(r),r),identity(3))
        g=F('9.80665')
        for sign in (1,-1):
            rr=[[F(1),F(0),F(0)],[F(0),F(3,5),-sign*F(4,5)],
                [F(0),sign*F(4,5),F(3,5)]]
            acc=[[F(0)],[sign*F(4,5)*g],[-F(3,5)*g]]
            self.assertEqual(matmul(transpose(rr),acc),[[F(0)],[F(0)],[-g]])
            self.assertEqual(matmul(transpose(rr),[[F(75)],[F(0)],[F(0)]]),
                             [[F(75)],[F(0)],[F(0)]])
        self.assertTrue(c['both_histories_nonzero_motion'])
        self.assertFalse(c['filter_divergence_claimed'])

    def test_service_is_all_phase_and_exceeds_contract(self):
        c=service_certificate()
        self.assertGreater(F(c['actual_innovation_service_lower']),3)
        self.assertGreater(F(c['true_heading_axis_service_lower']),1)
        self.assertGreater(F(c['all_phase_prefix_determinant_lower']),0)
        self.assertGreater(F(c['first_partial_cell_determinant_lower']),0)
        p=[[F(v) for v in row] for row in c['post_mag_covariance_upper']]
        self.assertTrue(is_psd(add(p,correct(predict(p,F('.04')),F('.8')**2/75**2),F(-1))))

    def test_batch_root_information_equals_actual_innovation_loss(self):
        # Independent dense, exact sequential-vs-batch check. Root observation
        # information survives interleaved prediction; endpoint covariance is
        # deliberately not substituted for root information.
        p0=[[F(1,20),F(1,100)],[F(1,100),F(1,10)]]
        p=[row[:] for row in p0]; m=identity(2); loss=[[F(0)]*2 for _ in range(2)]
        spacing=F(1,5); r=F(1,7); n=4
        obs=[]; unconditional=[]; cov=p0
        # Observation covariance entries from the unconditioned Markov model.
        for j in range(n):
            cov=predict(cov,spacing); unconditional.append(cov)
            obs.append([F(1),(j+1)*spacing])
        yy=[]
        for i in range(n):
            row=[]
            for j in range(n):
                small,big=min(i,j),max(i,j)
                cross=matmul(unconditional[small],transpose(transition((big-small)*spacing)))
                row.append(cross[0][0]+(r if i==j else 0))
            yy.append(row)
        for _ in range(n):
            p=predict(p,spacing); m=matmul(transition(spacing),m)
            s=p[0][0]+r; row=m[0][:]
            loss=add(loss,[[x*y/s for y in row] for x in row])
            gain=[p[i][0]/s for i in range(2)]
            m=[[m[i][j]-gain[i]*row[j] for j in range(2)] for i in range(2)]
            p=correct(p,r)
        batch=matmul(matmul(transpose(obs),inverse(yy)),obs)
        self.assertEqual(loss,batch)

if __name__=='__main__':
    unittest.main()
