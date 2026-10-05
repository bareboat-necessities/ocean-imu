import json,unittest
from pathlib import Path
from tools.stability.ou3_theorem.planar_service_verify import verify

ROOT=Path(__file__).resolve().parents[2]/'reports/results/ou3_stability'

class PlanarServiceVerifyTests(unittest.TestCase):
    def data(self):
        return [json.loads((ROOT/name).read_text()) for name in (
            'planar-service-stream-diagnostic.json','planar-service-operation-audit.json',
            'planar-service-frechet-diagnostic.json')]

    def test_committed_finite_evidence(self):
        result=verify(*self.data())
        self.assertTrue(result['finite_reproduction_pass'])
        self.assertFalse(result['all_time_magnetic_service_verified'])

    def test_rejects_missing_control(self):
        a,b,c=self.data();a.pop('untapped_complete_tail_samples_bitwise_equal')
        with self.assertRaises(ValueError):verify(a,b,c)

    def test_rejects_promotion(self):
        for index in range(3):
            rows=self.data();rows[index]['all_time_magnetic_service_verified']=True
            with self.assertRaises(ValueError):verify(*rows)

    def test_rejects_unlinked_stream(self):
        a,b,c=self.data();c['stream_sha256']='0'*64
        with self.assertRaises(ValueError):verify(a,b,c)

    def test_tangent_free_reports_cannot_bypass_source_binding(self):
        a,b,c=self.data();a['record_counts'].pop('9',None)
        a['driver_sha256']='0'*64
        with self.assertRaises(ValueError):verify(a,b,c)

    def test_rejects_nan_gain_and_missing_reset(self):
        a,b,c=self.data();c['parity_blocks']['odd']['relative_Frobenius_gain']=float('nan')
        with self.assertRaises(ValueError):verify(a,b,c)
        a,b,c=self.data();b['counts']['8']-=1
        with self.assertRaises(ValueError):verify(a,b,c)

if __name__=='__main__':unittest.main()
