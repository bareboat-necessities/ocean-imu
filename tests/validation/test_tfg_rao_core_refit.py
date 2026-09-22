"""Selection, sealed-holdout and failed-case retention contracts for RAO refit."""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import tfg_rao_core_refit as core
import tfg_rao_refit_analysis as analysis
import tfg_rao_refit_replay as replay


def fixture_rows():
    rows=[]
    for config in (core.arm('baseline'),core.arm('candidate',tau=1.025)):
        for record in replay.RECORDS:
            for seed in core.SEEDS['refine']:
                rows.append(dict(config=config['name'],env=config['env'],family='TFG',
                    input=record.filename,seed=seed,exit_code=0,violations=[],error=None,
                    metrics=dict.fromkeys(analysis.METRICS,1.0)))
    return rows


class SelectionTests(unittest.TestCase):
    def test_joint_refinement_is_four_dimensional(self):
        plan=json.loads((ROOT/'tools/tfg_rao_refit_plan.json').read_text())
        configs=core.configs_for('refine',plan)
        self.assertEqual(len(configs),82)
        points=[core.coeffs(c['env']) for c in configs[1:]]
        self.assertEqual(len(set(points)),81)
        self.assertEqual(sum(p[2]!=p[3] for p in points),54)
        self.assertEqual(len({p[3] for p in points if p[:3]==points[0][:3]}),3)
        sets=[set(core.SEEDS[p]) for p in ('screen','refine','holdout')]
        self.assertFalse(sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2])

    def test_gate_penalty_is_paired_not_net(self):
        rows=fixture_rows()
        base=[r for r in rows if r['config']=='baseline']
        cand=[r for r in rows if r['config']=='candidate']
        base[0]['violations']=[{'message':'ERROR: old','ratio':1.1}]
        cand[1]['violations']=[{'message':'ERROR: new','ratio':1.1}]
        row=next(r for r in core.rank(rows) if r['config']=='candidate')
        self.assertAlmostEqual(row['penalty'],.025)
        other=next(r for r in analysis.analyse(rows)['ranking'] if r['config']=='candidate')
        self.assertAlmostEqual(other['score'],row['score'])

    def test_missing_nonfinite_and_duplicates_do_not_disappear(self):
        rows=fixture_rows()
        rows[-1]['metrics']['pitch_rms_deg']=None
        self.assertFalse(next(r for r in core.rank(rows) if r['config']=='candidate')['eligible'])
        result=analysis.analyse(rows)
        self.assertEqual(len(result['ranking']),2)
        self.assertFalse(next(r for r in result['ranking'] if r['config']=='candidate')['admissible'])
        with self.assertRaises(ValueError): analysis.analyse(rows+[rows[0]])
        with self.assertRaises(ValueError): core.cases(rows+[rows[0]])
        with self.assertRaises(ValueError): analysis.analyse(rows[1:])

    def test_balanced_and_pooled_are_distinct(self):
        rows=fixture_rows()
        high=replay.RECORDS[-1].filename
        for row in rows:
            truth=1000. if row['input']==high else 1.
            ratio=(1.5 if row['input']==high else .9) if row['config']=='candidate' else 1.
            row['metrics']={m:truth*ratio for m in analysis.METRICS}
        result=next(r for r in analysis.analyse(rows)['ranking'] if r['config']=='candidate')
        self.assertLess(result['raw_score'],0)
        self.assertGreater(result['metrics']['disp_3d_rms_m']['pooled_change_pct'],0)
        self.assertAlmostEqual(result['raw_score'],(7*math.log(.9)+math.log(1.5))/8)

    def test_freeze_is_immutable_and_rejects_holdout_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            source=directory/'selection.json'
            rows=fixture_rows()
            source.write_text(json.dumps(rows))
            argv=['refit','--phase','freeze','--output-dir',str(directory),
                  '--selection-runs',str(source),'--candidate-config','candidate']
            with patch.object(sys,'argv',argv): core.main()
            frozen=(directory/'candidate.json').read_bytes()
            with patch.object(sys,'argv',argv),self.assertRaises(ValueError): core.main()
            self.assertEqual((directory/'candidate.json').read_bytes(),frozen)
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            source=directory/'selection.json'
            rows[0]['seed']=core.SEEDS['holdout'][0]
            source.write_text(json.dumps(rows))
            argv=['refit','--phase','freeze','--output-dir',str(directory),'--selection-runs',str(source)]
            with patch.object(sys,'argv',argv),self.assertRaises(ValueError): core.main()
            self.assertFalse((directory/'candidate.json').exists())


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.binary=self.root/'sim'
        self.records=[replay.RECORDS[0]]
        (self.root/self.records[0].filename).write_text('fixture input')
        self.configs=[{'name':'baseline','env':{}}]
        self.provenance={'input_sha256':{self.records[0].filename:'fixture-hash'}}

    def tearDown(self):
        self.temp.cleanup()

    def fake(self,code=0,nonfinite=False):
        metrics=dict.fromkeys(replay.METRICS,1.)
        metrics.update(samples=180000,window_s=900,family='TFG')
        if nonfinite: metrics['pitch_rms_deg']='nan'
        text='VALIDATION_METRICS '+' '.join(f'{k}={v}' for k,v in metrics.items())
        self.binary.write_text('#!/usr/bin/env python3\nimport os\nprint('+repr(text)+')\n'
            'print("QUALITY_GATE: PASS='+('1' if code==0 else '0')+'")\n'
            'print("PARENT_SEED="+str(os.environ.get("W3D_INIT_SEED")))\nraise SystemExit('+str(code)+')\n')
        self.binary.chmod(0o755)

    def execute(self,out):
        with patch.dict(replay.FAMILIES,{'TFG':self.binary}), \
             patch.object(replay,'RECORDS',self.records), \
             patch.object(replay,'input_provenance',return_value=self.provenance), \
             patch.object(replay,'source_commit',return_value='fixture'):
            return replay.execute(self.configs,'TFG',[0],['default'],self.root/out,1,self.root/'cache')

    def test_resume_checks_payload_and_clears_inherited_seeds(self):
        self.fake()
        with patch.dict('os.environ',{'W3D_INIT_SEED':'unwanted'}): first=self.execute('a')
        self.assertIsNone(first[0]['error'])
        self.assertTrue(self.execute('b')[0]['reused'])
        self.assertIn('PARENT_SEED=None',next((self.root/'a').glob('*.log')).read_text())
        cache=next((self.root/'cache').glob('*.json'))
        envelope=json.loads(cache.read_text())
        envelope['payload']['seconds']=999
        cache.write_text(json.dumps(envelope))
        with self.assertRaises(ValueError): self.execute('c')

    def test_fatal_and_nonfinite_cases_are_retained(self):
        for code,nonfinite in [(3,False),(0,True)]:
            with self.subTest(code=code,nonfinite=nonfinite):
                self.fake(code,nonfinite)
                rows=self.execute(f'failure-{code}')
                self.assertIsNotNone(rows[0]['error'])
                summary=json.loads((self.root/f'failure-{code}'/'summary.json').read_text())
                self.assertFalse(summary[0]['eligible'])
                self.assertEqual(len(rows),1)

    def test_gate_and_unexposed_overrides_are_rejected(self):
        with self.assertRaises(ValueError):
            replay.validate_configs([{'name':'b','env':{'TFG_GATE':'2'}}],b'TFG_GATE\0')
        with self.assertRaises(ValueError):
            replay.validate_configs([{'name':'b','env':{'TFG_TAU_COEFF':'2'}}],b'')


if __name__=='__main__':
    unittest.main()
