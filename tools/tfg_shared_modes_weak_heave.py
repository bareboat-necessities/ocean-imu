#!/usr/bin/env python3
"""Targeted diagnosis of the primary study's reduced-heave failure.

Exploratory, triggered by the registered primary result. Both geometries use
known coupling, real-only coupling, unchecked estimation and coupling-off
references. The operational estimator is not retuned on these outcomes.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,json,subprocess,tempfile
import numpy as np
import pandas as pd
from tfg_shared_modes_study import ROOT,OUT,BINARY,SCENARIOS,clean_env,sha,ci,save_trace
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--jobs',type=int,default=2)
ap.add_argument('--analysis-only',action='store_true',help='Reanalyse completed immutable simulation caches after an analysis-only edit.')
args=ap.parse_args()
filters=('A41@C','T41@C','CSM','TCSM','CSM:off','TCSM:off','CSM:noconf','TCSM:noconf',
 'CSM/real','TCSM/real','CSM:noconf/ext-true','TCSM:noconf/ext-true','TCSM/noreset')
seeds=list(range(1001,1013));cache=OUT/'weak_cache';cache.mkdir(exist_ok=True)
protocol=dict(binary_sha256=sha(BINARY),runner_sha256=sha(__file__),filters=filters,seeds=seeds,
 scenario=SCENARIOS['weak_heave']+'+T1200',role='exploratory failure diagnosis; no estimator retuning; primary multiplicity family unchanged',
 heading_diagnostic='RMS heading above 30 deg is reported separately from the inherited position/nonfinite divergence flag')
p=OUT/'weak_protocol.json'
if p.exists():
 old=json.loads(p.read_text());new=json.loads(json.dumps(protocol))
 if args.analysis_only:
  if {k:v for k,v in old.items() if k!='runner_sha256'}!={k:v for k,v in new.items() if k!='runner_sha256'}:
   raise RuntimeError('Simulation protocol changed; analysis-only reuse is invalid')
 elif old!=new:raise RuntimeError('Weak-heave protocol changed')
else:p.write_text(json.dumps(protocol,indent=2)+'\n')
def run(seed):
 target=cache/f'{seed}.json'
 if target.exists():return json.loads(target.read_text())
 if args.analysis_only:raise RuntimeError('Missing simulation cache in analysis-only mode')
 with tempfile.TemporaryDirectory(prefix='tfg-weak-heave-') as td:
  env=clean_env();env['CSM_TRACE']=str(Path(td)/'confidence_')
  cmd=[str(BINARY),'--header','--extended','--coupling','--level','full','--scenario',protocol['scenario'],
   '--seed',str(seed),'--score-start','300','--filters',','.join(filters)]
  tr=Path(td)/'trace.csv'
  if seed==1001:cmd+=['--timeseries',str(tr),'--decimate','20']
  r=subprocess.run(cmd,env=env,text=True,capture_output=True)
  if r.returncode:raise RuntimeError(r.stderr[-2000:])
  import csv,io
  rows=list(csv.DictReader(io.StringIO(r.stdout)));diags=[]
  for f in filters:
   path=Path(td)/('confidence_'+f.replace('/','_').replace(':','_')+'.csv')
   if path.exists():
    d=pd.read_csv(path);w=d[d.t>=300]
    diags.append(dict(filter=f,seed=seed,mean_confidence=float(w.confidence.mean()),
      min_confidence=float(w.confidence.min()),max_confidence=float(w.confidence.max()),shrinks=int(d.shrinks.iloc[-1]),
      mean_coherence_xz=float(w.coh_xz.mean()),mean_nis_ratio=float(w.nis_ratio.mean())))
    if seed==1001 and f in ('CSM','TCSM','CSM/real','TCSM/real'):
     save_trace(d,OUT/'traces'/('weak_confidence_'+f.replace('/','_')+'.csv.gz'))
  if seed==1001:
   d=pd.read_csv(tr);cols=['t']+[c for c in d if c.startswith(('dyaw_','droll_','dpitch_')) and c.split('_',1)[1] in ('A41@C','T41@C','CSM','TCSM','CSM/real','TCSM/real')]
   save_trace(d[cols],OUT/'traces'/'weak_attitude.csv.gz')
  record=dict(seed=seed,rows=rows,diagnostics=diags);target.write_text(json.dumps(record)+'\n');return record
records=[]
with ThreadPoolExecutor(max_workers=args.jobs) as pool:
 fs=[pool.submit(run,s) for s in seeds]
 for i,f in enumerate(as_completed(fs),1):records.append(f.result());print(f'{i}/{len(fs)}',flush=True)
rows=[dict(row,seed=r['seed']) for r in records for row in r['rows']];df=pd.DataFrame(rows)
for c in df:
 if c not in ('level','scenario','filter'):df[c]=pd.to_numeric(df[c],errors='coerce')
df=df.sort_values(['seed','filter']);df.to_csv(OUT/'weak_runs.csv',index=False,float_format='%.9g')
diag=pd.DataFrame([d for r in records for d in r['diagnostics']]).sort_values(['seed','filter']);diag.to_csv(OUT/'weak_confidence.csv',index=False,float_format='%.9g')
summary=df.groupby('filter').agg(n=('seed','size'),heave_m=('rms_z','mean'),horizontal_m=('rms_h','mean'),
 heading_deg=('yaw_rms_deg','mean'),heading_over30=('yaw_rms_deg',lambda x:int((x>30).sum())),
 tilt_deg=('tilt_rms_deg','mean'),nees9=('nees9_phys','mean'),nis=('nis','mean'),flagged=('failed','sum'))
summary=summary.join(diag.groupby('filter')[['mean_confidence','shrinks']].mean())
summary.to_csv(OUT/'weak_summary.csv',float_format='%.9g');paired=[]
for a,b in [('TCSM','CSM'),('TCSM:noconf/ext-true','CSM:noconf/ext-true'),('TCSM/real','CSM/real'),
 ('TCSM:off','CSM:off'),('TCSM/noreset','TCSM'),('CSM/real','CSM'),('TCSM/real','TCSM')]:
 aa=df[df['filter']==a].set_index('seed');bb=df[df['filter']==b].set_index('seed')
 for metric in ('rms_z','rms_h','yaw_rms_deg','tilt_rms_deg','nees9_phys'):
  n,d,lo,hi,pv=ci(aa[metric]-bb[metric]);paired.append(dict(filter=a,baseline=b,metric=metric,n=n,delta=d,ci_low=lo,ci_high=hi,p=pv))
pd.DataFrame(paired).to_csv(OUT/'weak_paired.csv',index=False,float_format='%.9g')
(OUT/'weak_analysis.json').write_text(json.dumps(dict(analysis_sha256=sha(__file__),execution_protocol_sha256=sha(p),analysis_only=args.analysis_only),indent=2)+'\n')
print(summary.to_string());print('Weak-heave diagnosis complete',flush=True)
