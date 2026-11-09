#!/usr/bin/env python3
"""Paired full-rate transition windows and explicit recovery/censoring criteria.

Uses the unchanged benchmark executable. Position CSV resolution is 10 um;
these diagnostics concern mm-and-larger transients, not sub-resolution gains.
Primary scores remain the benchmark's unrounded internal scores.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,json,subprocess,tempfile
import numpy as np
import pandas as pd
from tfg_shared_modes_study import ROOT,OUT,BINARY,SCENARIOS,PRIMARY,clean_env,sha,ci
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--jobs',type=int,default=2);args=ap.parse_args()
cases=('roundtrip','direction_step');seeds=list(range(1001,1025));cache=OUT/'transient_cache';cache.mkdir(exist_ok=True)
protocol=dict(binary_sha256=sha(BINARY),runner_sha256=sha(__file__),seeds=seeds,cases=cases,
 filters=PRIMARY,sample_rate_hz=200,position_export_resolution_m=1e-5,
 criterion='first 20 consecutive 1-Hz evaluations of trailing 10-s RMS below 1.25 times the post-transition steady-window RMS; time measured from transition end; otherwise censored',
 events={'roundtrip':{'up':[410,780,600,780],'down':[810,1200,1000,1200]},'direction_step':{'direction':[601.5,1200,900,1200]}},
 role='secondary paired transient diagnostic; primary inference unchanged')
p=OUT/'transient_protocol.json'
if p.exists() and json.loads(p.read_text())!=json.loads(json.dumps(protocol)):raise RuntimeError('Transient protocol changed')
p.write_text(json.dumps(protocol,indent=2)+'\n')
def rms(x):return float(np.sqrt(np.mean(np.asarray(x)**2)))
def recovery(t,sq,start,stop,ref_start,ref_stop):
 reference=np.mean(sq[(t>=ref_start)&(t<ref_stop)]);threshold=1.25*np.sqrt(reference)
 # Exact 200 Hz box sums, then 1 Hz evaluations; never score drawing traces.
 power=np.r_[0.,np.cumsum(np.nan_to_num(sq))];i=np.arange(2000,len(t)+1,200)
 rt=t[i-1];rr=np.sqrt(np.maximum(0.,(power[i]-power[i-2000])/2000))
 valid=(rt>=start)&(rt<stop);rt=rt[valid];under=rr[valid]<=threshold
 for j in range(max(0,len(rt)-19)):
  if under[j:j+20].all():return float(rt[j]-start),0,float(threshold)
 return float('nan'),1,float(threshold)
def run(task):
 case,seed=task;target=cache/f'{case}_{seed}.json'
 if target.exists():return json.loads(target.read_text())
 with tempfile.TemporaryDirectory(prefix='tfg-transients-') as td:
  trace=Path(td)/'trace.csv';cmd=[str(BINARY),'--header','--extended','--level','full',
   '--scenario',SCENARIOS[case]+'+T1200','--seed',str(seed),'--score-start','300','--filters',','.join(PRIMARY),
   '--timeseries',str(trace),'--decimate','1']
  proc=subprocess.run(cmd,env=clean_env(),capture_output=True,text=True)
  if proc.returncode:raise RuntimeError(proc.stderr[-2000:])
  cols=lambda c:c=='t' or c.startswith(('x_','y_','z_'))
  d=pd.read_csv(trace,usecols=cols);t=d.t.to_numpy();rows=[]
  for f in PRIMARY:
   ep=np.stack([d[a+'_'+f]-d[a+'_true'] for a in 'xyz'],axis=1);ez=ep[:,2]**2;eh=np.sum(ep[:,:2]**2,axis=1)
   for event,(end,stop,r0,r1) in protocol['events'][case].items():
    lo=end-(20 if case=='roundtrip' else 3);window=(t>=lo)&(t<end+60)
    zr,zc,zt=recovery(t,ez,end,stop,r0,r1);hr,hc,ht=recovery(t,eh,end,stop,r0,r1)
    rows.append(dict(case=case,seed=seed,filter=f,event=event,rms_z=float(np.sqrt(np.mean(ez[window]))),
      rms_h=float(np.sqrt(np.mean(eh[window]))),peak_z=float(np.sqrt(np.max(ez[window]))),peak_h=float(np.sqrt(np.max(eh[window]))),
      recovery_z_s=zr,recovery_h_s=hr,censored_z=zc,censored_h=hc,threshold_z_m=zt,threshold_h_m=ht))
  target.write_text(json.dumps(rows)+'\n');return rows
rows=[]
with ThreadPoolExecutor(max_workers=args.jobs) as pool:
 futures=[pool.submit(run,(c,s)) for c in cases for s in seeds]
 for i,f in enumerate(as_completed(futures),1):rows.extend(f.result());print(f'{i}/{len(futures)}',flush=True)
df=pd.DataFrame(rows).sort_values(['case','seed','filter','event']);df.to_csv(OUT/'transient_runs.csv',index=False,float_format='%.9g')
paired=[]
for (case,event),g in df.groupby(['case','event']):
 b=g[g['filter']=='CSM'].set_index('seed')
 for f in PRIMARY:
  if f=='CSM':continue
  a=g[g['filter']==f].set_index('seed')
  for metric in ['rms_z','rms_h','peak_z','peak_h','recovery_z_s','recovery_h_s']:
   n,d,lo,hi,pv=ci(a[metric]-b[metric]);paired.append(dict(case=case,event=event,filter=f,metric=metric,n=n,delta=d,ci_low=lo,ci_high=hi,p=pv,censored_z=int(a.censored_z.sum()),censored_h=int(a.censored_h.sum())))
pd.DataFrame(paired).to_csv(OUT/'transient_paired.csv',index=False,float_format='%.9g');print('Transient diagnostics complete',flush=True)
