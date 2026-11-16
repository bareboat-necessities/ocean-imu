"""Publication never accepts a wrong case, short study, changed file or header."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import shared_modes_article_evidence as study


class SharedModesArticleEvidenceTests(unittest.TestCase):
    def fixture(self, root):
        results = root/'results'
        results.mkdir()
        (root/'src').mkdir()
        (root/'src/model.h').write_text('shipping header')
        subprocess.run(['git', 'init', '-q', str(root)], check=True)
        subprocess.run(['git', '-C', str(root), 'add', 'src'], check=True)
        subprocess.run(['git', '-C', str(root), '-c', 'user.name=Test', '-c', 'user.email=test@example.org', 'commit', '-qm', 'source'], check=True)
        revision = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'],text=True).strip()
        with (results/'scores.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['key','score_samples']);writer.writeheader()
            writer.writerows({'key':str(i),'score_samples':180000} for i in range(1038))
        protocol = dict(dataset=dict(input_sha256=study.sim_dataset.REFERENCE_CSV_SHA256),
                        seed_triplets=list(range(10)),source_hashes={},source_revision=revision)
        (results/'protocol.json').write_text(json.dumps(protocol))
        (results/'case-receipt.json').write_text(json.dumps(dict(case=study.sim_dataset.CASE,status='completed',mode='full',source_commit=revision,inputs=protocol['dataset'])))
        (results/'engine_summary.csv').write_text('retained engine scores')
        summary=dict(replays=1038,files_sha256={p.name:study.sha(p) for p in (results/'scores.csv',results/'protocol.json',results/'engine_summary.csv')})
        (results/'publication_summary.json').write_text(json.dumps(summary))
        return results

    def test_valid_complete_study(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);study.validate(self.fixture(root),root)

    def test_wrong_case_or_failed_receipt_rejected(self):
        for field,value in [('case','invalid'),('status','failed'),('source_commit','0'*40)]:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as t:
                root=Path(t);r=self.fixture(root)
                p=r/'case-receipt.json';d=json.loads(p.read_text());d[field]=value;p.write_text(json.dumps(d))
                with self.assertRaises(ValueError):study.validate(r,root)

    def test_short_window_rejected_even_if_product_hash_updated(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);r=self.fixture(root)
            p=r/'scores.csv';p.write_text(p.read_text().replace('180000','12000',1))
            p=r/'publication_summary.json';s=json.loads(p.read_text());s['files_sha256']['scores.csv']=study.sha(r/'scores.csv');p.write_text(json.dumps(s))
            with self.assertRaisesRegex(ValueError,'shortened'):study.validate(r,root)

    def test_uncommitted_runtime_header_change_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);r=self.fixture(root);(root/'src/model.h').write_text('changed')
            with self.assertRaisesRegex(ValueError,'replay sources changed'):study.validate(r,root)

    def test_changed_input_or_product_rejected(self):
        for name in ('protocol.json','scores.csv'):
            with self.subTest(name=name),tempfile.TemporaryDirectory() as t:
                root=Path(t);r=self.fixture(root)
                if name=='protocol.json':
                    p=r/name;d=json.loads(p.read_text());d['dataset']['input_sha256']={};p.write_text(json.dumps(d))
                else:
                    p=r/name;p.write_text(p.read_text().replace('0,180000','x,180000'))
                with self.assertRaises(ValueError):study.validate(r,root)
