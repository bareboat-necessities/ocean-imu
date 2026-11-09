#!/usr/bin/env python3
"""Diagnostic reset-transport supplement, reusing the registered data runner.

This includes the existing T41R0 research variant as well as the new shared
variant. It does not change the primary 2x2 estimands or their multiplicity set.
"""
import argparse,json,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
from tfg_shared_modes_study import ROOT,OUT,BINARY,AB_SCENARIOS,execute_one,ci,sha
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--jobs',type=int,default=1);args=ap.parse_args()
filters=('T41@C','T41R0@C','TCSM','TCSM/noreset');seeds=list(range(1001,1013))
protocol=dict(binary_sha256=sha(BINARY),runner_sha256=sha(ROOT/'tools/tfg_shared_modes_study.py'),
 source_sha256=sha(__file__),filters=filters,seeds=seeds,cases=AB_SCENARIOS,role='exploratory reset-transport diagnostic; paired and independent of the primary inference family')
p=OUT/'reset_protocol.json'
if p.exists() and json.loads(p.read_text())!=json.loads(json.dumps(protocol)):raise RuntimeError('Reset supplement protocol changed')
p.write_text(json.dumps(protocol,indent=2)+'\n')
tasks=[('reset',n,s,filters) for n in AB_SCENARIOS for s in seeds]
todo=[t for t in tasks if not (OUT/'cache'/f'{t[0]}_{t[1]}_{t[2]}.json').exists()]
with ThreadPoolExecutor(max_workers=args.jobs) as pool:
 futures=[pool.submit(execute_one,t,OUT) for t in todo]
 for i,f in enumerate(as_completed(futures),1):r=f.result();print(f'{i}/{len(todo)} {r["key"]}',flush=True)
rows=[]
for task in tasks:
 rec=json.loads((OUT/'cache'/f'{task[0]}_{task[1]}_{task[2]}.json').read_text())
 for row in rec['rows']:rows.append(dict(case=rec['case'],**row))
df=pd.DataFrame(rows)
for c in df:
 if c not in ('level','scenario','filter','case'):df[c]=pd.to_numeric(df[c],errors='coerce')
df.to_csv(OUT/'reset_runs.csv',index=False,float_format='%.9g')
rows=[]
for case,g in df.groupby('case'):
 for a,b,geometry in [('TCSM/noreset','TCSM','Shared'),('T41R0@C','T41@C','Independent')]:
  aa=g[g['filter']==a].set_index('seed');bb=g[g['filter']==b].set_index('seed')
  for m in ['rms_z','rms_h','rms_3d','tilt_rms_deg','yaw_rms_deg','nees9_phys','nis']:
   n,d,lo,hi,pv=ci(100*(aa[m]/bb[m]-1));rows.append(dict(case=case,geometry=geometry,metric=m,n=n,no_transport_pct=d,ci_low=lo,ci_high=hi,p=pv,failures=int(aa.failed.sum())))
pd.DataFrame(rows).to_csv(OUT/'reset_comparison.csv',index=False,float_format='%.9g')
print('Reset supplement complete',flush=True)
