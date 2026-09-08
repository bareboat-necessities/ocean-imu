#!/usr/bin/env python3
"""Evaluate explicit filter settings on the pinned RAO records; never alter gates."""
import argparse,hashlib,json,os,re,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from model_mismatch_ablation import FAMILIES,RECORDS,source_commit
from ou_validation import parse_validation_metrics
from sim_dataset import input_provenance
ROOT=Path(__file__).resolve().parents[1]
def main():
 producer_bytes=Path(__file__).read_bytes()
 revision=source_commit()
 source_diff=subprocess.check_output(['git','diff','HEAD','--','src','tools','tests'],cwd=ROOT)
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--configs',type=Path,required=True);p.add_argument('--family',default='OU-III',choices=FAMILIES)
 p.add_argument('--records',default='0,3,4,7');p.add_argument('--seeds',default='default');p.add_argument('--jobs',type=int,default=4)
 p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 configs=json.loads(a.configs.read_text()); records=[RECORDS[int(i)] for i in a.records.split(',')]
 binary=FAMILIES[a.family];a.output_dir.mkdir(parents=True,exist_ok=True)
 binary_bytes=binary.read_bytes()
 paths={r.filename:(binary.parent/r.filename).resolve() for r in records}
 prov=input_provenance(paths.values())
 # Performance controls only. Reject any attempted test/gate or sensor-data override.
 for c in configs:
  for k in c['env']:
   if (not k.startswith(('SF_','OU_II_','OU_III_','TFG_')) and k not in ('W3D_AW_COV_SYNC','OU_SIGMA_STILL_DECAY_SEC','OU_SIGMA_VAR_K_PERIODS','OU_SIGMA_VAR_HORIZON_MIN_S','OU_SIGMA_VAR_HORIZON_MAX_S')) or any(x in k for x in ('GATE','LIMIT','SEED')):raise ValueError(k)
   # A family may not expose another family's knob. Require its exact C-string
   # in the executable so an ignored variable cannot masquerade as an ablation.
   # This is a necessary check; paired controls still verify the applied effect.
   if k.encode()+b'\0' not in binary_bytes:raise ValueError(f'{k} is not exposed by {binary}')
 tasks=[(c,r,s) for c in configs for r in records for s in a.seeds.split(',')]
 def run(task):
  c,r,s=task;key=f"{c['name']}_{r.spectrum}_{r.hs_m}_{s}";dst=a.output_dir/(key+'.json')
  env=dict(os.environ,W3D_WRITE_TIMESERIES='0',W3D_COLLECT_ALL_GATES='1',W3D_VALIDATION_WINDOW_SEC='900')
  if s!='default':env.update(W3D_INIT_SEED=s,W3D_IMU_SEED=s)
  env.update(c['env']);start=time.monotonic()
  q=subprocess.run([str(binary),'--input',str(paths[r.filename])],cwd=binary.parent,env=env,capture_output=True,text=True)
  (a.output_dir/(key+'.log')).write_text(q.stdout+q.stderr)
  if q.returncode not in (0,1) or (q.returncode==1 and 'QUALITY_GATE: PASS=0' not in q.stdout):raise RuntimeError(key+' '+str(q.returncode))
  violations=[]
  for line in (q.stdout+q.stderr).splitlines():
   if line.startswith('ERROR:'):
    m=re.search(r'\(([0-9.eE+-]+)[^>]*>\s*([0-9.eE+-]+)',line)
    violations.append(dict(message=line,ratio=float(m[1])/float(m[2]) if m else None))
  row=dict(config=c['name'],env=c['env'],family=a.family,input=r.filename,seed=s,exit_code=q.returncode,seconds=time.monotonic()-start,violations=violations,metrics=parse_validation_metrics(q.stdout))
  dst.write_text(json.dumps(row,indent=2)+'\n');return row
 rows=[]
 with ThreadPoolExecutor(a.jobs) as pool:
  for f in as_completed([pool.submit(run,t) for t in tasks]):
   row=f.result();rows.append(row);print(f"{len(rows)}/{len(tasks)} {row['config']} {row['input']} seed={row['seed']} gates={row['exit_code']} yaw={row['metrics'].get('yaw_rms_deg')}",flush=True)
 (a.output_dir/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
 if binary.read_bytes()!=binary_bytes:raise RuntimeError('Simulator changed during the study')
 (a.output_dir/'manifest.json').write_text(json.dumps(dict(source_commit=revision,source_dirty=bool(source_diff),source_diff_sha256=hashlib.sha256(source_diff).hexdigest(),simulation_provenance=prov,configs=configs,records=a.records,seeds=a.seeds,producer_sha256=hashlib.sha256(producer_bytes).hexdigest(),binary_sha256=hashlib.sha256(binary_bytes).hexdigest(),protocol='Full records, final 900 s; unchanged executable quality gates; seeded paired sensor and initialization draws'),indent=2)+'\n')
 summary=[]
 for c in configs:
  rr=[r for r in rows if r['config']==c['name']];vv=[v['ratio'] for r in rr for v in r['violations'] if v['ratio'] is not None]
  summary.append(dict(config=c['name'],failed_records=sum(r['exit_code']!=0 for r in rr),violations=sum(len(r['violations']) for r in rr),worst_ratio=max(vv,default=1),**{k:sum(r['metrics'][k] for r in rr)/len(rr) for k in ['yaw_rms_deg','roll_rms_deg','pitch_rms_deg','disp_3d_rms_m','accel_bias_3d_rms_mps2']}))
 (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(sorted(summary,key=lambda r:(r['violations'],r['worst_ratio'])),indent=2))
if __name__=='__main__':main()
