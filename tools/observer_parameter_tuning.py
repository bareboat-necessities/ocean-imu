#!/usr/bin/env python3
"""Paired NLO/PII parameter replay with immutable executable quality gates."""
import argparse,hashlib,json,os,re,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from model_mismatch_ablation import RECORDS,source_commit
from sim_dataset import input_provenance
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=['nlo','pii'],required=True)
 p.add_argument('--configs',type=Path,required=True);p.add_argument('--seeds',default='default');p.add_argument('--jobs',type=int,default=2);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 folder='nlo' if a.family=='nlo' else 'pii_observer';binary=ROOT/'tests'/folder/('nlo-sim' if a.family=='nlo' else 'pii_observer-adaptive')
 configs=json.loads(a.configs.read_text());a.output_dir.mkdir(parents=True,exist_ok=True)
 paths={r.filename:(binary.parent/r.filename).resolve() for r in RECORDS};prov=input_provenance(paths.values());binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
 for c in configs:
  for k in c['env']:
   if k not in ['PII_R_SCALE','PII_TAU_A_SCALE','PII_KP_SCALE','PII_KI_SCALE','NLO_THETA','NLO_THETA_GAIN','NLO_THETA_TAU','NLO_COARSE_TAU','NLO_K1I','NLO_KII']:raise ValueError(k)
 def run(task):
  c,seed=task;key=c['name']+'_'+seed
  with tempfile.TemporaryDirectory(prefix='observer-tuning-') as tmp:
   for r in RECORDS:
    (Path(tmp)/r.filename).symlink_to(paths[r.filename])
    suffix='_tvg_nlo_nomag_nognss.csv' if a.family=='nlo' else '_nonkalman_fusion.csv'
    (Path(tmp)/r.filename.replace('wave_data_','w3d_').replace('.csv',suffix)).symlink_to('/dev/null')
   env=dict(os.environ,W3D_COLLECT_ALL_GATES='1');env.update(c['env'])
   if seed!='default':env.update(W3D_IMU_SEED=seed,W3D_INIT_SEED=seed)
   q=subprocess.run([str(binary)],cwd=tmp,env=env,capture_output=True,text=True)
   text=q.stdout+q.stderr;(a.output_dir/(key+'.log')).write_text(text)
   if q.returncode not in (0,1) or (q.returncode==1 and 'QUALITY_GATE: PASS=0' not in text):raise RuntimeError(key+' '+str(q.returncode))
   z=[float(v) for v in re.findall(r'Z RMS(?: raw)? \(%Hs\): ([0-9.eE+-]+)',text)]
   angles=[list(map(float,v)) for v in re.findall(r'Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw(?:\(free/ignored\))?=([0-9.eE+-]+)',text)]
   if len(z)!=8 or len(angles)!=8:raise RuntimeError('Missing records: '+key)
   row=dict(config=c['name'],env=c['env'],seed=seed,exit_code=q.returncode,z_pct_hs=z,angles_deg=angles,failures=[l for l in text.splitlines() if l.startswith('ERROR:')])
   (a.output_dir/(key+'.json')).write_text(json.dumps(row,indent=2)+'\n');print(key,len(row['failures']),flush=True);return row
 with ThreadPoolExecutor(a.jobs) as pool:rows=list(pool.map(run,[(c,s) for c in configs for s in a.seeds.split(',')]))
 if hashlib.sha256(binary.read_bytes()).hexdigest()!=binary_hash:raise RuntimeError('Simulator changed during experiment')
 (a.output_dir/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
 (a.output_dir/'manifest.json').write_text(json.dumps(dict(source_commit=source_commit(),simulation_provenance=prov,binary_sha256=binary_hash,producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),configs=configs,seeds=a.seeds,protocol='Full eight records, final 900 s; paired IMU and initialization seeds; no quality gate changes'),indent=2)+'\n')
if __name__=='__main__':main()
