#!/usr/bin/env python3
"""Paired 2x2 TFG/shared-CARMA experiment. Full-rate scores, bounded drawing data.

Run: python3 tools/tfg_shared_modes_study.py run --jobs 5
Then: python3 tools/tfg_shared_modes_study.py analyse
No oracle channel enters an operational filter. Named ext-true/flip variants
are explicitly diagnostic. Cache identity includes executable and source hashes.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,os,subprocess,sys,tempfile,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'plots'))
from vector_export import envelope_indices
OUT=ROOT/'reports/results/tfg_shared_modes'
BINARY=ROOT/'tests/kalman_carma_shared_modes/kalman_carma_shared_modes-bench'
BASE='jonswap_fp0.29_hs0.8_g3.3'
PRIMARY=('A41@C','T41@C','CSM','TCSM')
SCENARIOS={
 'jonswap':BASE,
 'pm_long':'pm_fp0.16_hs1.5',
 'regular':'regular_f0.25_h1',
 'short_low':'jonswap_fp0.40_hs0.27_g3.3',
 'high_long':'jonswap_fp0.12_hs4_g3.3',
 'broad_direction':BASE+'+sp1',
 'long_crested':BASE+'+sp10000',
 'heading90':BASE+'+rot90',
 'head_seas':BASE+'+rel180',
 'swell':'swell',
 'mixed':'crossing',
 'weak_heave':BASE+'+hscale0.15',
 'calm':'calm',
 'low_energy':'small_amplitude',
 'engine800':BASE+'+engine800',
 'engine1600':BASE+'+engine1600',
 'engine2400':BASE+'+engine2400',
 'engine3200':BASE+'+engine3200',
 'colored_noise':BASE+'+cnoise4',
 'roundtrip':BASE+'+start0.5+lev400w20x4+lev800w20x0.5',
 'direction_step':'dir_change+rapid_direction',
 'changing_period':'changing_fp',
 'maneuver':BASE+'+maneuver+turns',
 'rapid_rotation':'large_rotation+turns',
 'mag_outage':BASE+'+maggap300+gbdrift',
 'bias_temperature':'bias_drift+gbdrift+temp',
 'jitter':BASE+'+jit',
}
ABLATIONS=('TCSM','CSM','TCSM:off','TCSM/real','CSM/real','TCSM/noreset',
            'TCSM:noconf/ext-true','TCSM:noconf/ext-true/real',
            'TCSM/ext-flip','TCSM:noconf/ext-flip')
AB_SCENARIOS=('long_crested','mixed','direction_step','roundtrip','rapid_rotation','mag_outage')
METRICS=('rms_z','rms_h','rms_3d','rmsv_x','rmsv_y','rmsv_z','rmsv_h',
         'roll_rms_deg','pitch_rms_deg','tilt_rms_deg','yaw_rms_deg','bg_rms_dph','ba_rms_mg',
         'nis','nis_in95','nees_pos3','nees_pv6','att_nees','nees9_phys','nees9_native',
         'cov95_x','cov95_y','cov95_z','su_rms_z','su_peak_z','su_tilt_rms_deg','su_yaw_rms_deg',
         'conv_z','conv_x','conv_y','yaw_gap_rms_deg','yaw_post_rms_deg','yaw_recover_s',
         'turn_peak_deg','turn_recover_s')
SOURCES=('src/kalman_carma/Kalman3D_Wave_TFG_CARMA3.h',
         'src/kalman_carma/Kalman3D_Wave_CARMA3.h',
         'src/kalman_carma_shared_modes/CarmaSharedModesMotion.h',
         'src/kalman_carma_shared_modes/SharedModesEstimator.h',
         'src/kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h',
         'tests/kalman_carma/carma3_bench.cpp','tests/kalman_carma/CarmaBenchTruth.h',
         'tests/kalman_carma_shared_modes/CsmBenchAdapter.h','tools/tfg_shared_modes_study.py')
def sha(p):
 """Hash exact file bytes for experiment source and binary provenance."""
 return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean_env():
 """Remove inherited filter overrides before running controlled simulations."""
 return {k:v for k,v in os.environ.items() if not k.startswith(('W3D_','CSM_','OU_','SF_','TFG_'))}
def save_trace(frame,path):
 """Preserve joint extrema in bounded traces shared by SVG and PGF exports."""
 # A single shared extrema selection bounds both SVG and PGF geometry.
 idx=envelope_indices(tuple(frame[c].to_numpy() for c in frame),max_points=2000)
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as z:
  z.write(frame.iloc[idx].to_csv(index=False,float_format='%.8g').encode())
def execute_one(task,out):
 """Run one paired scenario/seed cell and cache scores plus selected traces."""
 stage,name,seed,filters=task;sc=SCENARIOS[name]+'+T1200'
 key=f'{stage}_{name}_{seed}'
 trace=seed==1001 and name in ('jonswap','roundtrip','direction_step','mag_outage','rapid_rotation','engine2400') and stage=='primary'
 with tempfile.TemporaryDirectory(prefix='tfg-shared-') as td:
  tr=Path(td)/'trace.csv'
  cmd=[str(BINARY),'--header','--extended','--coupling','--level','full','--scenario',sc,
       '--seed',str(seed),'--score-start','300','--filters',','.join(filters)]
  if trace:cmd+=['--timeseries',str(tr),'--decimate','1']
  start=time.monotonic();p=subprocess.run(cmd,env=clean_env(),text=True,capture_output=True)
  if p.returncode:raise RuntimeError(f'{key}: {p.returncode}\n{p.stdout[-1000:]}\n{p.stderr[-2000:]}')
  rows=list(csv.DictReader(io.StringIO(p.stdout)))
  if {r['filter'] for r in rows}!=set(filters):raise RuntimeError(f'{key}: missing filter rows')
  if trace:
   df=pd.read_csv(tr)
   # Store families independently: selecting extrema jointly over >100 channels
   # can undersample narrow transients. Every family has a 2000-point bound.
   for family,prefix in {'position':('x_','y_','z_'),'velocity':('vx_','vy_','vz_'),
                          'attitude':('droll_','dpitch_','dyaw_','roll_','pitch_','yaw_'),
                          'gyro_bias':('bgx_','bgy_','bgz_'),'accel_bias':('bax_','bay_','baz_'),
                          'uncertainty':('ex_','ey_','ez_','sx_','sy_','sz_')}.items():
    cols=['t']+[c for c in df if c.startswith(prefix)]
    save_trace(df[cols],out/'traces'/f'{name}_{family}.csv.gz')
    if name=='jonswap':save_trace(df.loc[df.t>=1170,cols],out/'traces'/f'{name}_tail_{family}.csv.gz')
   # Full-rate error windows for degradation and recovery, not drawing data.
   segments={'pre':(300,390),'rising':(390,440),'high':(500,780),'falling':(790,840),'recovered':(950,1200)} if name=='roundtrip' else {}
   if name=='direction_step':segments={'pre':(400,590),'turn':(598.5,620),'post60':(603,663),'recovered':(900,1200)}
   segment_rows=[]
   for label,(lo,hi) in segments.items():
    win=df[(df.t>=lo)&(df.t<hi)]
    for f in filters:
     e=np.stack([win[f'{a}_{f}']-win[f'{a}_true'] for a in 'xyz'],axis=1)
     segment_rows.append(dict(scenario=name,seed=seed,segment=label,filter=f,
           rms_z=float(np.sqrt(np.mean(e[:,2]**2))),rms_h=float(np.sqrt(np.mean(np.sum(e[:,:2]**2,axis=1))))))
   if segment_rows:pd.DataFrame(segment_rows).to_csv(out/'traces'/f'{name}_segments.csv',index=False)
  record=dict(key=key,stage=stage,case=name,seed=seed,command=cmd[:],elapsed_s=time.monotonic()-start,rows=rows,stderr=p.stderr[-3000:])
  path=out/'cache'/f'{key}.json';path.parent.mkdir(parents=True,exist_ok=True)
  path.write_text(json.dumps(record,separators=(',',':'))+'\n')
  return record

def run(args):
 """Register a fresh study protocol and execute all missing paired cells."""
 out=args.results;out.mkdir(parents=True,exist_ok=True)
 protocol=dict(baseline_revision='1fe4711b489ad92f3eef8ab560da904af741df9d',
   binary_sha256=sha(BINARY),source_hashes={p:sha(ROOT/p) for p in SOURCES},
   sample_rate_hz=200,duration_s=1200,score_window_s=[300,1200],
   primary_filters=PRIMARY,primary_seeds=list(range(1001,1001+args.seeds)),
   ablation_filters=ABLATIONS,ablation_seeds=list(range(1001,1001+args.ablation_seeds)),
   scenarios=SCENARIOS,ablation_scenarios=AB_SCENARIOS,
   wave_excitation='Evidence for every CARMA cell; @C explicitly selects it on independent filters',
   design='paired 2x2; identical truth/IMU/magnetometer/initialization per seed and scenario; moving startup in every non-calm case',
   primary_endpoints=['rms_z','rms_h'],practical_improvement_percent=5.,
   inference='paired seed-level Student t intervals; marginal CIs, Holm-adjusted p across all primary scenario/endpoints; interaction is paired difference-in-differences',
   caveats='synthetic RAOs and sensor model, not sea trials; no ESP32 hardware timing')
 pp=out/'protocol.json'
 if pp.exists() and json.loads(pp.read_text())!=json.loads(json.dumps(protocol)):raise RuntimeError('Protocol changed. Choose a fresh results directory; do not mix cached experiments.')
 pp.write_text(json.dumps(protocol,indent=2)+'\n')
 tasks=[('primary',n,s,PRIMARY) for n in SCENARIOS for s in protocol['primary_seeds']]
 tasks += [('ablation',n,s,ABLATIONS) for n in AB_SCENARIOS for s in protocol['ablation_seeds']]
 todo=[t for t in tasks if not (out/'cache'/f'{t[0]}_{t[1]}_{t[2]}.json').exists()]
 print(f'{len(tasks)-len(todo)}/{len(tasks)} cached; {len(todo)} remaining',flush=True)
 with ThreadPoolExecutor(max_workers=args.jobs) as pool:
  futures={pool.submit(execute_one,t,out):t for t in todo}
  for i,future in enumerate(as_completed(futures),1):
   rec=future.result();bad=[r['filter'] for r in rec['rows'] if r['failed']!='0']
   print(f'{i}/{len(todo)} {rec["key"]} {rec["elapsed_s"]:.1f}s failed={bad}',flush=True)
 collect(out)

def collect(out):
 """Collect cached scores or retain the published compact scores on clean checkouts."""
 rows=[]
 for p in sorted((out/'cache').glob('*.json')):
  r=json.loads(p.read_text())
  for row in r['rows']:rows.append(dict(stage=r['stage'],case=r['case'],**row))
 if not rows:
  # Publication bundles retain compact scores, not restart caches. Reanalysis
  # must not truncate those scores to an empty CSV on a clean checkout.
  for name in ('runs.csv','runs.csv.gz'):
   if (out/name).exists():return pd.read_csv(out/name)
  raise RuntimeError('No experiment caches or published per-run scores in '+str(out))
 df=pd.DataFrame(rows)
 for c in df:
  if c not in ('level','scenario','filter','stage','case'):df[c]=pd.to_numeric(df[c],errors='coerce')
 # Runtime measurements are host-dependent; retain scores and runtime but never
 # compare timings in regression gates. Small per-run summaries are reproducible.
 df.to_csv(out/'runs.csv',index=False,float_format='%.9g')
 return df

def ci(x):
 """Compute a paired Student interval and zero-mean test over finite seed differences."""
 x=np.asarray(x,float);x=x[np.isfinite(x)];n=len(x)
 if n<2:return n,float('nan'),float('nan'),float('nan'),float('nan')
 mean=float(np.mean(x));se=float(stats.sem(x));half=float(stats.t.ppf(.975,n-1)*se)
 p=float(stats.ttest_1samp(x,0).pvalue) if se>0 else (0. if mean else 1.)
 return n,mean,mean-half,mean+half,p

def analyse(args):
 """Derive paired comparisons, multiplicity corrections and condition-balanced effects."""
 out=args.results;df=collect(out);primary=df[df.stage=='primary'];summ=[];paired=[];inter=[]
 comparisons=(('TFG geometry, shared','TCSM','CSM'),('Shared prior, TFG','TCSM','T41@C'),
              ('TFG geometry, independent','T41@C','A41@C'),('Shared prior, QMEKF','CSM','A41@C'))
 for (case,f),g in primary.groupby(['case','filter']):
  row=dict(case=case,filter=f,n=len(g),failures=int(g.failed.sum()))
  for m in METRICS:row[m]=g[m].mean()
  summ.append(row)
 for case,g in primary.groupby('case'):
  for label,a,b in comparisons:
   aa=g[g['filter']==a].set_index('seed');bb=g[g['filter']==b].set_index('seed');aa,bb=aa.align(bb,join='inner',axis=0)
   for m in METRICS:
    n,d,lo,hi,p=ci(aa[m]-bb[m]);nr,r,rl,rh,rp=ci(100*(aa[m]/bb[m]-1))
    paired.append(dict(case=case,comparison=label,metric=m,n=n,delta=d,ci_low=lo,ci_high=hi,p=p,
                       percent=r,pct_ci_low=rl,pct_ci_high=rh,pct_p=rp,
                       failed_a=int(aa.failed.sum()),failed_b=int(bb.failed.sum())))
  v=g.pivot(index='seed',columns='filter',values=['rms_z','rms_h','rms_3d','tilt_rms_deg','yaw_rms_deg'])
  for m in v.columns.levels[0]:
   x=(v[m]['TCSM']-v[m]['CSM'])-(v[m]['T41@C']-v[m]['A41@C'])
   n,d,lo,hi,p=ci(x);inter.append(dict(case=case,metric=m,n=n,interaction=d,ci_low=lo,ci_high=hi,p=p))
 pd.DataFrame(summ).to_csv(out/'summary.csv',index=False,float_format='%.9g')
 pairs=pd.DataFrame(paired);pairs['p_holm']=np.nan
 mask=(pairs.comparison=='TFG geometry, shared')&pairs.metric.isin(['rms_z','rms_h'])
 selected=pairs.loc[mask,'p'].sort_values();adjust=np.maximum.accumulate(np.minimum(1.,selected.to_numpy()*(len(selected)-np.arange(len(selected)))))
 pairs.loc[selected.index,'p_holm']=adjust
 pairs.to_csv(out/'paired.csv',index=False,float_format='%.9g');pd.DataFrame(inter).to_csv(out/'interactions.csv',index=False,float_format='%.9g')
 # Aggregate by seed first, giving each scenario equal weight; seeds, not
 # individual time samples, are the independent resampling units.
 aggregate=[]
 for m in ('rms_z','rms_h','rms_3d','tilt_rms_deg','yaw_rms_deg','bg_rms_dph','ba_rms_mg'):
  v=primary.pivot(index=['case','seed'],columns='filter',values=m)
  for label,a,b in comparisons:
   lg=np.log(v[a]/v[b]).groupby('seed').mean();n,d,lo,hi,p=ci(lg)
   aggregate.append(dict(metric=m,comparison=label,n=n,geomean_pct=100*np.expm1(d),ci_low_pct=100*np.expm1(lo),ci_high_pct=100*np.expm1(hi),p=p))
 pd.DataFrame(aggregate).to_csv(out/'aggregate.csv',index=False,float_format='%.9g')
 ab=df[df.stage=='ablation'];rows=[]
 for case,g in ab.groupby('case'):
  base=g[g['filter']=='TCSM'].set_index('seed')
  for f,x in g.groupby('filter'):
   if f=='TCSM':continue
   x=x.set_index('seed');x,b=x.align(base,join='inner',axis=0)
   for m in ('rms_z','rms_h','tilt_rms_deg','yaw_rms_deg','nees9_phys','nis'):
    n,d,lo,hi,p=ci(100*(x[m]/b[m]-1));rows.append(dict(case=case,filter=f,metric=m,n=n,percent=d,ci_low=lo,ci_high=hi,p=p,failures=int(x.failed.sum())))
 pd.DataFrame(rows).to_csv(out/'ablations.csv',index=False,float_format='%.9g')
 print(pd.DataFrame(aggregate).to_string(index=False));print('Failures',df.groupby(['stage','filter']).failed.sum().to_dict())
 # Stable report inputs; source data/cache can be regenerated and is not committed.
 print(f'Wrote {len(df)} per-run rows, {len(pairs)} paired comparisons.')

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['run','analyse','collect']);ap.add_argument('--results',type=Path,default=OUT)
 ap.add_argument('--jobs',type=int,default=5);ap.add_argument('--seeds',type=int,default=24);ap.add_argument('--ablation-seeds',type=int,default=12);args=ap.parse_args()
 {'run':run,'analyse':analyse,'collect':lambda a:collect(a.results)}[args.stage](args)
