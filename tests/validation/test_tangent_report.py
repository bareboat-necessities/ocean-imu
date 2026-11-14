"""Reject mixed evidence and represent unavailable statistical tests honestly."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from tangent_shared_modes_inputs import verify_raw_inputs

try:
    from tangent_shared_modes_report import ljung_box, trace
except ModuleNotFoundError:
    ljung_box = trace = None


class RawInputBindingTests(unittest.TestCase):
    def fixture(self, work):
        name = 'wave_17_native'
        contents = {'scores': 'scenario,seed,filter,start,end,sensor_hash,failed\n'
                    'wave,17,TFG,150,900,paired,0\n',
                    'trace': 't,filter,white_x\n899,TFG,0.25\n', 'model': '{}\n'}
        receipt = {'command': ['producer', '--scenario', 'wave', '--seed', '17',
                              '--output', '/original/computer/'+name]}
        for kind, content in contents.items():
            path = work/(name+'_'+kind+('.json' if kind == 'model' else '.csv'))
            path.write_text(content)
            receipt[kind+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        scores = work/'scores.csv'
        # The aggregate may have canonical LF endings and a different row order.
        scores.write_text('scenario,seed,filter,start,end,sensor_hash,failed,noise_profile\n'
                          'wave,17,TFG,150,900,paired,0,native\n')
        return {'commands': [receipt], 'seeds': [17]}, scores

    def test_relocated_verified_inputs_and_custom_seed(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            manifest, scores = self.fixture(work)
            self.assertEqual(verify_raw_inputs(work, manifest, scores), 1)
            if trace is not None:
                self.assertEqual(trace(work, seed=manifest['seeds'][0]).white_x.iloc[0], .25)

    def test_mismatched_raw_artifacts_are_rejected(self):
        for suffix in ('_trace.csv', '_scores.csv', '_model.json'):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                manifest, scores = self.fixture(work)
                path = work/('wave_17_native'+suffix)
                path.write_text(path.read_text()+'\n')
                with self.assertRaisesRegex(ValueError, 'raw input differs'):
                    verify_raw_inputs(work, manifest, scores)

    def test_published_scores_must_match_receipt_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            manifest, scores = self.fixture(work)
            scores.write_text(scores.read_text().replace('paired,0', 'paired,1'))
            with self.assertRaisesRegex(ValueError, 'published scores differ'):
                verify_raw_inputs(work, manifest, scores)


@unittest.skipIf(ljung_box is None, 'statistical reporting dependencies unavailable')
class InnovationWhitenessTests(unittest.TestCase):
    def test_short_nonfinite_and_constant_innovations_are_unavailable(self):
        for values in ([], [float('nan')]*30, list(range(20)), [2.]*30, [.1]*30):
            with self.subTest(values=values):
                result = ljung_box(values)
                self.assertIsNone(result['statistic'])
                self.assertIsNone(result['p'])
                json.dumps(result, allow_nan=False)

    def test_nonconstant_innovations_have_a_finite_statistic(self):
        result = ljung_box([(-1.)**i for i in range(60)])
        self.assertEqual(result['samples'], 60)
        self.assertGreater(result['statistic'], 0)
        self.assertLessEqual(result['p'], 1)
        json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
