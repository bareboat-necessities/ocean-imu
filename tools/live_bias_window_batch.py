#!/usr/bin/env python3
"""Retain seed-level window diagnostics and prior-sensitivity controls."""
import concurrent.futures,json,argparse
from pathlib import Path
from live_bias_window import run
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/results/live_bias_consistency'
def main():
 p=argparse.ArgumentParser();p.add_argument('--jobs',type=int,default=2);a=p.parse_args()
 tasks=[]
 for seed in range(1001,1007):
  for f in ('qmekf','tfg'):
   for motion in ('wave','rocking'):
    prefix=OUT/'runs/screen'/f'{seed}_{f}_{motion}_exact_fixed_both_none_none_1200'
    for ps,bs,label in ((1,.0452949761388,'native_prior'),(4,.0452949761388,'motion4'),(1,.1358849284164,'bias3')):
     tasks.append((prefix,ps,bs,label,seed,f,motion))
 def work(t):
  prefix,ps,bs,label,seed,f,motion=t
  x=run(prefix,ps,20,bs);x.update(seed=seed,filter=f,motion=motion,case=label)
  x['prefix']=str(prefix.relative_to(ROOT))
  return x
 rows=[]
 with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
  for i,x in enumerate(pool.map(work,tasks)):
   rows.append(x);(OUT/'window_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(i+1,len(tasks),x['status'],x['relinearized']['ba_h_mm_s2'],flush=True)
if __name__=='__main__':main()
