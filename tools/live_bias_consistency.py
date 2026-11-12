#!/usr/bin/env python3
"""Paired current-main controls, independent freezes, injections, long tails.

All runs reuse exactly the same raw noise/physical history per seed and motion.
No cache reuse after binary/protocol changes. Full-rate scores; bounded plots.
"""
import argparse,concurrent.futures,hashlib,json,re,subprocess,time,shutil
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/results/live_bias_consistency'
BIN=ROOT/'tests/live_bias_consistency/live-bias-followup'
def main():
    p=argparse.ArgumentParser();p.add_argument('--jobs',type=int,default=4);p.add_argument('--stage',choices=['screen','long','recovery','noise','baseline','startup','latency','calibration','combined','late_recovery'],default='screen');p.add_argument('--seeds',type=int,nargs='+',default=list(range(1001,1007)));p.add_argument('--label');p.add_argument('--binary',type=Path,default=BIN);p.add_argument('--noise',choices=['native','gyro','mag','all','prior','late_prior','all_prior','schmidt','schmidt_prior'],default='native');a=p.parse_args()
    binary=a.binary
    stage_label=a.label or (a.stage if a.noise=='native' else a.stage+'_'+a.noise)
    records=[]
    for seed in a.seeds:
      for f in ('qmekf','tfg'):
        if a.stage=='late_recovery':
          for motion in ('wave','rocking'):
            for injection in ('none','bias','tilt'):
              records.append((seed,f,motion,'exact_fixed','both_start','none',injection,7200))
          for arm in ('learned_fixed','correct_reference'):
            records.append((seed,f,'wave',arm,'both_start','none','none',7200))
        elif a.stage=='combined':
          for arm in ('exact_fixed','learned_fixed'):
            for noise in ('native','schmidt','schmidt_prior'):
              records.append((seed,f,'wave',arm,'both','none','none',7200,noise))
        elif a.stage=='calibration':
          for arm in ('exact_fixed','learned_fixed'):
            for noise in ('prior','late_prior','all_prior','schmidt','schmidt_prior'):
              records.append((seed,f,'wave',arm,'both','none','none',3600,noise))
        elif a.stage=='latency':
          for lag in (0,.005,.02):records.append((seed,f,'wave','exact_fixed','both','none','none',1200,'native',lag))
        elif a.stage=='startup':
          for timing in ('endpoint','both'):records.append((seed,f,'wave','learned_fixed',timing,'none','none',1200))
        elif a.stage=='noise':
          for motion in ('wave','rocking'):
            for timing in ('endpoint','both'):
              for noise in ('native','gyro','mag','all'):
                records.append((seed,f,motion,'exact_fixed',timing,'none','none',1200,noise))
        elif a.stage=='screen':
          for motion in ('wave','rocking'):
            for timing in ('endpoint','core','frontend','both'):
              for arm in ('exact_fixed','oracle'):
                records.append((seed,f,motion,arm,timing,'none','none',1200))
            for freeze in ('sea','coupling','both'):
              records.append((seed,f,motion,'exact_fixed','both',freeze,'none',1200))
          for timing in ('endpoint','both'):records.append((seed,f,'wave','learned_fixed',timing,'none','none',1200))
        elif a.stage in ('long','baseline'):
          for timing in (('endpoint',) if a.stage=='baseline' else ('endpoint','both')):
            for arm in ('exact_fixed','oracle'):
              records.append((seed,f,'wave',arm,timing,'none','none',7200))
          records.append((seed,f,'rocking','exact_fixed','endpoint' if a.stage=='baseline' else 'both','none','none',3600))
        else:
          for injection in ('bias','tilt'):
            for motion in ('wave','rocking'):
              records.append((seed,f,motion,'exact_fixed','both','none',injection,3600))
          records.append((seed,f,'wave','correct_reference','both','none','none',3600))
    if a.stage not in ('noise','latency','calibration','combined'):records=[tuple(r)+(a.noise,) for r in records]
    if a.stage!='latency':records=[tuple(r)+(0.,) for r in records]
    runroot=OUT/'runs'/stage_label;runroot.mkdir(parents=True,exist_ok=True)
    source_binary=binary.resolve();snapshot=runroot/'estimator'
    shutil.copy2(source_binary,snapshot);binary=snapshot.resolve()
    binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    def run(row):
      seed,f,motion,arm,timing,freeze,inj,duration,noise,lag=row
      key=f'{seed}_{f}_{motion}_{arm}_{timing}_{freeze}_{inj}_{duration}';key+=('_'+noise if noise!='native' else '');key+=('_lag'+str(lag) if lag else '');prefix=runroot/key
      cmd=[str(binary),'--noise',noise,'--seed',str(seed),'--filter',f,'--motion',motion,'--arm',arm,'--timing',timing,'--freeze',freeze,'--injection',inj,'--duration',str(duration),'--sigma','.001','--output',str(prefix)]
      if lag:cmd+=['--mag-lag',str(lag)]
      if a.stage in ('startup','combined'):cmd+=['--change','0']
      if a.stage=='late_recovery':cmd+=['--change','3600']
      if a.stage=='screen' and timing=='both' and arm=='exact_fixed' and freeze=='none':cmd+=['--window','1']
      start=time.monotonic();r=subprocess.run(cmd,text=True,capture_output=True);prefix.with_suffix('.log').write_text(r.stdout+r.stderr)
      if r.returncode:raise RuntimeError(key+': '+r.stderr)
      parsed=dict(re.findall(r'(\w+)=([^\s]+)',r.stdout));scores=pd.read_csv(str(prefix)+'_scores.csv');diag=pd.read_csv(str(prefix)+'_diag.csv')
      meta=dict(seed=seed,filter=f,motion=motion,arm=arm,timing=timing,freeze=freeze,injection=inj,duration=duration,noise=noise,lag=lag,stage=stage_label,seconds=time.monotonic()-start,**parsed)
      results=[]
      for w in scores.to_dict('records'):
        dd=diag[(diag.t>=w['start'])&(diag.t<w['end'])]
        results.append(dict(**meta,**w,coverage=float((dd.h_nees<=5.99146).mean()),h_nees=float(dd.h_nees.mean()),h_sigma=float(dd.h_sigma.mean()),nuisance_inflation=float(dd.nuisance_inflation.mean())))
      return results
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
      for i,res in enumerate(pool.map(run,records)):
        rows+=res
        pd.DataFrame(rows).to_csv(OUT/f'{stage_label}_runs.csv',index=False,float_format='%.10g')
        print(f'{a.stage} {i+1}/{len(records)}',flush=True)
    df=pd.DataFrame(rows)
    for _,group in df.groupby(['seed','motion','duration','lag']):
      assert group.original_sensor_hash.nunique()==1,'unpaired sensors'
    manifest=dict(baseline='9018fa4b57d9bd10983ae4e77ff95ce509551e17',binary_sha256=binary_hash,stage=a.stage,runs=len(records),seeds=a.seeds,paired=True,commands=' '.join(__import__('sys').argv))
    (OUT/f'{stage_label}_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
