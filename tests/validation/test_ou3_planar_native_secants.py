import tempfile
import unittest
from pathlib import Path
import numpy as np
from tools.stability.ou3_theorem.planar_native_secants import difference, svec, read, analyze

class NativeSecantTests(unittest.TestCase):
    def test_left_quaternion_chart_and_antipodal_representation(self):
        ref=np.zeros(466);ref[21:25]=[np.cos(.2),0,np.sin(.2),0]
        state=ref.copy();a=.001
        # q=Rx(a) q_ref, so left relative log is a e_x, not q_ref^-1 Rx q_ref.
        state[21:25]=[np.cos(a/2)*np.cos(.2),np.sin(a/2)*np.cos(.2),
                      np.cos(a/2)*np.sin(.2),np.sin(a/2)*np.sin(.2)]
        np.testing.assert_allclose(difference(state,ref)[:3],[a,0,0],atol=1e-15)
        state[21:25]*=-1
        np.testing.assert_allclose(difference(state,ref)[:3],[a,0,0],atol=1e-15)
    def test_symmetric_covariance_coordinates_preserve_frobenius_norm(self):
        rng=np.random.default_rng(653); A=rng.normal(size=(21,21));A=(A+A.T)/2
        self.assertAlmostEqual(np.linalg.norm(A,'fro'),np.linalg.norm(svec(A)),places=12)
    def test_partial_evidence_and_invalid_steps_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad';p.write_bytes(b'not a complete snapshot')
            with self.assertRaises(ValueError):read(p)
            for eps in (0,-.01,float('nan'),float('inf')):
                with self.assertRaises(ValueError):analyze(d,eps)
if __name__=='__main__':unittest.main()
