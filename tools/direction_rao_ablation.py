#!/usr/bin/env python3
"""End-to-end RAO direction/engine-guard ablation on the pinned vessel records."""
import argparse,csv,hashlib,json,os,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from model_mismatch_ablation import FAMILIES,RECORDS,source_commit
from ou_validation import parse_validation_metrics
from sim_dataset import input_provenance
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--jobs',type=int,default=3);p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/direction_rao_ablation');a=p.parse_args()
 paths={r.filename:(a.data_dir/r.filename).resolve() for r in RECORDS};prov=input_provenance(paths.values());a.output_dir.mkdir(parents=True,exist_ok=True)
 tasks=[(f,r,rpm,arm) for f in ['OU-II','OU-III'] for r in RECORDS for rpm in [0,2400] for arm in ['unconditioned','guard-only','guard-and-rao']]
 def run(task):
  family,r,rpm,arm=task;binary=FAMILIES[family]
  env=dict(os.environ,W3D_WRITE_TIMESERIES='0',W3D_VALIDATION_WINDOW_SEC='900',W3D_COLLECT_ALL_GATES='1',W3D_DIRECTION_RAO='vessel-rao-28ft' if arm=='guard-and-rao' else 'off')
  if arm=='unconditioned':env.update({('OU_III' if family=='OU-III' else 'OU_II')+suffix:'0' for suffix in ['_ACC_GUARD_HZ','_ACC_GUARD_RACC_GAIN']})
  if rpm:env.update(W3D_ENGINE_RPM=str(rpm),W3D_ENGINE_LEVEL_MPS2='0.6',W3D_ENGINE_BANDWIDTH_HZ='80')
  q=subprocess.run([str(binary),'--input',str(paths[r.filename])],cwd=binary.parent,env=env,capture_output=True,text=True)
  if q.returncode not in (0,1) or (q.returncode==1 and 'QUALITY_GATE: PASS=0' not in q.stdout):raise RuntimeError(str(task)+'\n'+q.stderr[-2000:])
  if rpm and 'ENGINE_VIBRATION' not in q.stdout:raise RuntimeError('Engine model did not engage')
  if ('ACC_GUARD' in q.stdout)!=(arm!='unconditioned'):raise RuntimeError('Incorrect guard arm')
  metrics=parse_validation_metrics(q.stdout)
  row=dict(metrics,family=family,input=r.filename,hs_m=r.hs_m,rpm=rpm,arm=arm,regression_exit_code=q.returncode)
  print(f'{family} {r.filename} rpm={rpm} {arm} axis={metrics.get("dir_axis_abs_error_deg")} unresolved={metrics.get("dir_sense_uncertain_pct")}',flush=True)
  return row
 with ThreadPoolExecutor(a.jobs) as pool:rows=list(pool.map(run,tasks))
 output=a.output_dir/'direction_runs.csv'
 with output.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (a.output_dir/'manifest.json').write_text(json.dumps(dict(source_commit=source_commit(),simulation_provenance=prov,protocol='Full eight records x two OU families x engine off/2400 RPM (0.6 m/s2 reference level, 80 Hz sensor bandwidth) x three conditioning arms; default draw; final 900 s; unchanged gates',producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),binary_sha256={k:hashlib.sha256(FAMILIES[k].read_bytes()).hexdigest() for k in ['OU-II','OU-III']},files={output.name:hashlib.sha256(output.read_bytes()).hexdigest()}),indent=2)+'\n')
if __name__=='__main__':main()
