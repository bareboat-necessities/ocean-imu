import struct,tempfile,unittest
from pathlib import Path
import numpy as np
from tools.stability.ou3_theorem.planar_service_stream import MAGIC,SIZES,records,expand,instrument,HEADER
from tools.stability.ou3_theorem.planar_parity import EVEN,ODD

class PlanarServiceStreamTests(unittest.TestCase):
    def read(self,content):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'stream';p.write_bytes(content);return list(records(p))
    def record(self,kind,sample,a):
        return struct.pack('<III',kind,sample,len(a))+np.asarray(a,dtype='<f4').tobytes()
    def test_every_operation_layout_is_accepted(self):
        data=MAGIC+b''.join(self.record(k,i+1,np.zeros(n)) for i,(k,n) in enumerate(SIZES.items()))
        self.assertEqual(len(self.read(data)),9)
    def test_bad_or_partial_evidence_fails_closed(self):
        bad=[b'badmagic',MAGIC+b'\x00',MAGIC+struct.pack('<III',1,1,2),
             MAGIC+self.record(1,1,np.zeros(SIZES[1]))[:-1],
             MAGIC+self.record(7,2,np.zeros(SIZES[7]))+self.record(7,1,np.zeros(SIZES[7])),
             MAGIC+self.record(7,1,np.full(SIZES[7],np.nan))]
        for data in bad:
            with self.subTest(size=len(data)),self.assertRaises(ValueError):self.read(data)
    def test_parity_expansion_preserves_cross_terms_within_blocks(self):
        a=np.arange(225,dtype=float);P=expand(a)
        np.testing.assert_array_equal(P[np.ix_(EVEN,EVEN)],a[:144].reshape(12,12))
        np.testing.assert_array_equal(P[np.ix_(ODD,ODD)],a[144:].reshape(9,9))
        np.testing.assert_array_equal(P[np.ix_(EVEN,ODD)],np.zeros((12,9)))
    def test_source_taps_bind_actual_pre_prediction_and_noise(self):
        s=instrument(HEADER.read_text())
        self.assertIn('planar_before_prediction(Pext.data())',s)
        self.assertIn('Pext.data(),R_S.data(),pseudo_update_elapsed_s_',s)
        self.assertIn('planar_after_sync(Pext.data(),aw_covariance_floor_target_.data(),planar_sync_pending)',s)
        self.assertIn('planar_after_reset(Pext.data())',s)
        self.assertIn('planar_mean_tangent(\"acc\"',s)
        self.assertIn('planar_mean_tangent(\"mag\"',s)
        with self.assertRaises(ValueError):instrument('wrong source')

if __name__=='__main__':unittest.main()
