#!/usr/bin/env python3
"""Reduced-heave dose response: the shared-mode excitation cap, before and after.

The heave RAO is scaled by k with horizontal motion unchanged, on the first 12
primary seeds. `/cap0` rows restore the legacy cap (the cross term scaled
down to applied coherence 0.999); unsuffixed shared rows are the operational
estimators, which keep the measured regression and enlarge the horizontal
prior at the cap (CarmaCouplingConfig::cap_coherence = 0.8).

Legacy mechanism, registered in sequence on the legacy estimator:
  Dose        shared heading error grows as the heave RAO shrinks, while
              real-only and independent-axis controls stay near nominal;
  Spreading   at k = 0.15 broad spreading weakens and long-crested seas
              strengthen the collapse;
  Cap         the known prior's applied coherence follows 0.334 / k^2 and is
              pinned at the cap for k < 0.58, where heave harm begins;
  Intervention  S_factor 2 moves that boundary to k = 0.29.
The cap policy was then changed: its target was chosen on seeds 1001-1003
of the k = 1, 0.5, 0.3 and 0.15 cells (0.7-0.8 equivalent, 0.95 and the
measured-coherence shrinkage worse). These cells compare it with the legacy
cap and both independent-axis filters on identical data. Known coupling is a
diagnostic input only; results are outside the primary multiplicity family.

Run:        python3 tools/tfg_shared_modes_heave_sweep.py run --jobs 3
Reanalyse:  python3 tools/tfg_shared_modes_heave_sweep.py analyse
"""
from __future__ import annotations
import argparse,csv,io,json,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from tfg_shared_modes_study import OUT,BINARY,BASE,clean_env,sha,ci

OPERATIONAL=('A41@C','T41@C','CSM','TCSM')
LEGACY=('CSM/cap0','TCSM/cap0','CSM/real/cap0','TCSM/real/cap0','CSM:noconf/ext-true/cap0')
FILTERS=OPERATIONAL+LEGACY+('CSM:noconf/ext-true',)
KNOWN_LEGACY,KNOWN_FIXED='CSM:noconf/ext-true/cap0','CSM:noconf/ext-true'
SEEDS=list(range(1001,1013))
# (cell, heave RAO factor, spreading exponent s; None keeps the default s=4)
CELLS=(('h1.00',1.,None),('h0.80',.8,None),('h0.70',.7,None),('h0.60',.6,None),('h0.50',.5,None),('h0.40',.4,None),('h0.30',.3,None),('h0.25',.25,None),
       ('h0.20',.2,None),('h0.15',.15,None),('h0.15_broad',.15,1),('h0.15_long',.15,10000))
DOSE=[c for c,_,s in CELLS if s is None]
CAP=.999
FIXED_CAP=.8
INTERVENTION=dict(s_factor=2.,cells=('h1.00','h0.50','h0.30','h0.25','h0.20'),
 filters=('A41@C/sf2','CSM/cap0/sf2','TCSM/cap0/sf2','CSM/real/cap0/sf2','TCSM/real/cap0/sf2','CSM:noconf/ext-true/cap0/sf2'))
HEADING_FAILURE_DEG=30.
REPRODUCE={'h1.00':('runs.csv.gz','jonswap'),'h0.15':('weak_runs.csv',None)}
SCORES=('rms_z','rms_h','yaw_rms_deg','yaw_mean_deg','tilt_rms_deg','ba_rms_mg','bg_rms_dph','nees9_phys','nis','failed',
        'cpl_coh_xz','cpl_ph_xz','cpl_coh_yz','cpl_ph_yz')

def scenario(scale,spread):
 """Compose the registered nominal sea with the heave and spreading modifiers."""
 s=BASE+('' if scale==1. else f'+hscale{scale:g}')+('' if spread is None else f'+sp{spread}')
 return s+'+T1200'

def protocol():
 """Execution identity: any change requires a fresh results directory."""
 return dict(binary_sha256=sha(BINARY),filters=FILTERS,seeds=SEEDS,score_start_s=300,
  cells={c:scenario(h,s) for c,h,s in CELLS},heading_failure_deg=HEADING_FAILURE_DEG,
  role='exploratory dose response; legacy-cap rows are paired controls; outside the primary multiplicity family',
  intervention=dict(INTERVENTION,prediction='legacy cap with S_factor 2 reaches the cap only for k < 0.29'),
  legacy_predictions={'dose':'shared heading error increases as heave RAO decreases; real-only and independent controls stay nominal',
   'spreading':'at k = 0.15 broad spreading weakens and long-crested seas strengthen the shared heading collapse',
   'cap':'the known prior is pinned at the 0.999 cap for k < 0.58; shared heave harm relative to real-only coupling begins only there'},
  fixed_policy=dict(cap_coherence=FIXED_CAP,development='seeds 1001-1003 of the k = 1, 0.5, 0.3 and 0.15 cells'))

def run_cell(stage,cell,scale,spread,seed,cache,filters):
 """Run one paired process (every filter sees the same truth and sensors)."""
 target=cache/(f'{cell}_{seed}.json' if stage=='dose' else f'{stage}_{cell}_{seed}.json')
 if target.exists():return json.loads(target.read_text())
 with tempfile.TemporaryDirectory(prefix='tfg-heave-sweep-') as td:
  env=clean_env();env['CSM_TRACE']=str(Path(td)/'confidence_')
  cmd=[str(BINARY),'--header','--extended','--coupling','--level','full','--scenario',scenario(scale,spread),
       '--seed',str(seed),'--score-start','300','--filters',','.join(filters)]
  p=subprocess.run(cmd,env=env,text=True,capture_output=True)
  if p.returncode:raise RuntimeError(f'{cell}/{seed}: {p.returncode}\n{p.stderr[-2000:]}')
  rows=list(csv.DictReader(io.StringIO(p.stdout)))
  if {r['filter'] for r in rows}!=set(filters):raise RuntimeError(f'{cell}/{seed}: missing filter rows')
  diags=[]
  for f in filters:
   path=Path(td)/('confidence_'+f.replace('/','_').replace(':','_')+'.csv')
   if path.exists():
    w=pd.read_csv(path);w=w[w.t>=300]
    diags.append(dict(filter=f,mean_confidence=float(w.confidence.mean()),min_confidence=float(w.confidence.min())))
 record=dict(stage=stage,cell=cell,seed=seed,command=cmd,rows=rows,diagnostics=diags)
 target.write_text(json.dumps(record,separators=(',',':'))+'\n')
 return record

def run(args):
 """Register the protocol and execute missing paired processes."""
 out=args.results;cache=out/'heave_sweep_cache';cache.mkdir(parents=True,exist_ok=True)
 path=out/'heave_sweep_protocol.json';new=json.loads(json.dumps(protocol()))
 if path.exists() and json.loads(path.read_text())!=new:
  raise RuntimeError('Heave-sweep protocol changed; use a fresh results directory')
 path.write_text(json.dumps(new,indent=2)+'\n')
 spec={c:(h,s) for c,h,s in CELLS}
 tasks=[('dose',c,h,s,seed,cache,FILTERS) for c,h,s in CELLS for seed in SEEDS]
 tasks+=[('sf2',c,*spec[c],seed,cache,INTERVENTION['filters']) for c in INTERVENTION['cells'] for seed in SEEDS]
 with ThreadPoolExecutor(max_workers=args.jobs) as pool:
  futures=[pool.submit(run_cell,*t) for t in tasks]
  for i,f in enumerate(as_completed(futures),1):
   r=f.result();print(f'{i}/{len(tasks)} {r["cell"]} {r["seed"]}',flush=True)
 analyse(args)

def load(out):
 """Collect the compact per-run scores from restart caches or the committed table."""
 recs=[json.loads(p.read_text()) for p in sorted((out/'heave_sweep_cache').glob('*.json'))]
 if not recs:return pd.read_csv(out/'heave_sweep_runs.csv')
 rows=[]
 for r in recs:
  conf={d['filter']:d for d in r['diagnostics']}
  for row in r['rows']:
   d=dict(stage=r.get('stage','dose'),cell=r['cell'],seed=r['seed'],filter=row['filter'],**{k:row[k] for k in SCORES})
   d['mean_confidence']=conf.get(row['filter'],{}).get('mean_confidence',np.nan)
   rows.append(d)
 df=pd.DataFrame(rows)
 for c in df:
  if c not in ('stage','cell','filter'):df[c]=pd.to_numeric(df[c],errors='coerce')
 return df.sort_values(['stage','cell','seed','filter']).reset_index(drop=True)

def reproduction(out,df):
 """Exact agreement of repeated registered seeds with the published per-run tables."""
 res={}
 for cell,(name,case) in REPRODUCE.items():
  ref=pd.read_csv(out/name)
  if case is not None:ref=ref[(ref.case==case)&(ref.stage=='primary')]
  ref=ref[ref.seed.isin(SEEDS)&ref['filter'].isin(FILTERS)].set_index(['seed','filter'])
  new=df[df.cell==cell].set_index(['seed','filter']).loc[ref.index]
  cols=['rms_z','rms_h','yaw_rms_deg','tilt_rms_deg','nees9_phys']
  diff=float(np.nanmax(np.abs(new[cols].to_numpy()-ref[cols].to_numpy())/np.maximum(np.abs(ref[cols].to_numpy()),1e-12)))
  res[cell]=dict(reference=name,rows=int(len(ref)),filters=sorted(set(ref.index.get_level_values(1))),max_relative_difference=diff)
 return res

def analyse(args):
 """Legacy mechanism, the fixed cap policy and its comparison with every filter."""
 out=args.results;runs=load(out);runs.to_csv(out/'heave_sweep_runs.csv',index=False,float_format='%.9g')
 df=runs[runs.stage=='dose'].drop(columns='stage');iv=runs[runs.stage=='sf2'].drop(columns='stage')
 if len(df)!=len(CELLS)*len(SEEDS)*len(FILTERS) or len(iv)!=len(INTERVENTION['cells'])*len(SEEDS)*len(INTERVENTION['filters']):
  raise RuntimeError('Incomplete heave sweep')
 summary=runs.groupby(['stage','cell','filter']).agg(n=('seed','size'),heave_mm=('rms_z',lambda x:1000*x.mean()),
  horizontal_mm=('rms_h',lambda x:1000*x.mean()),heading_deg=('yaw_rms_deg','mean'),heading_median_deg=('yaw_rms_deg','median'),
  heading_over30=('yaw_rms_deg',lambda x:int((x>HEADING_FAILURE_DEG).sum())),tilt_deg=('tilt_rms_deg','mean'),
  ba_mg=('ba_rms_mg','mean'),nees9=('nees9_phys','median'),flagged=('failed','sum'),
  applied_coherence_xz=('cpl_coh_xz','mean'),mean_confidence=('mean_confidence','mean'))
 summary.to_csv(out/'heave_sweep_summary.csv',float_format='%.9g')
 def pct(c,f,b,m):
  """Paired percent change of filter f against b in one cell."""
  v=df[df.cell==c].pivot(index='seed',columns='filter',values=m);n,d,lo,hi,_=ci(100*(v[f]/v[b]-1))
  return dict(percent=d,ci_low=lo,ci_high=hi)
 paired=[]
 for cell,g in df.groupby('cell'):
  for metric in ('yaw_rms_deg','rms_z','rms_h'):
   v=g.pivot(index='seed',columns='filter',values=metric)
   for a,b in (('CSM','A41@C'),('TCSM','A41@C'),('TCSM','T41@C'),('TCSM','CSM'),('CSM','CSM/cap0'),('TCSM','TCSM/cap0'),
               ('CSM/cap0','A41@C'),('TCSM/cap0','A41@C'),('TCSM/cap0','CSM/cap0'),
               ('CSM/cap0','CSM/real/cap0'),('TCSM/cap0','TCSM/real/cap0')):
    n,d,lo,hi,p=ci(v[a]-v[b]);paired.append(dict(cell=cell,filter=a,baseline=b,metric=metric,n=n,delta=d,ci_low=lo,ci_high=hi,p=p))
 for f in FILTERS:
  for other in ('h0.15_broad','h0.15_long'):
   v=df[df['filter']==f].pivot(index='seed',columns='cell',values='yaw_rms_deg')
   n,d,lo,hi,p=ci(v[other]-v['h0.15']);paired.append(dict(cell=other+'-h0.15',filter=f,baseline=f,metric='yaw_rms_deg',n=n,delta=d,ci_low=lo,ci_high=hi,p=p))
 pr=pd.DataFrame(paired);pr.to_csv(out/'heave_sweep_paired.csv',index=False,float_format='%.9g')
 s=summary.reset_index();s=s[s.stage=='dose']
 def mean(cell,f,col='heading_deg'):return float(s[(s.cell==cell)&(s['filter']==f)][col].iloc[0])
 onset={}
 for f in ('CSM/cap0','TCSM/cap0','CSM','TCSM'):
  hit=[c for c in DOSE if pr[(pr.cell==c)&(pr['filter']==f)&(pr.baseline=='A41@C')&(pr.metric=='yaw_rms_deg')].ci_low.iloc[0]>1.]
  onset[f]=hit[0] if hit else None
 mc=df.assign(mc=df.cpl_coh_xz+df.cpl_coh_yz).groupby(['cell','filter']).mc.mean()
 cap={c:dict(known_legacy=float(mc[(c,KNOWN_LEGACY)]),known_fixed=float(mc[(c,KNOWN_FIXED)]),
            estimated_legacy=float(mc[(c,'CSM/cap0')]),estimated_fixed=float(mc[(c,'CSM')]),
            known_at_cap=bool(mc[(c,KNOWN_LEGACY)]>=CAP-5e-4)) for c in DOSE}
 harm={}
 for f,b in (('CSM/cap0','CSM/real/cap0'),('TCSM/cap0','TCSM/real/cap0')):
  harm[f]={c:pct(c,f,b,'rms_z') for c in DOSE}
  hit=[c for c in DOSE if harm[f][c]['ci_low']>0];harm[f]['first_harm_cell']=hit[0] if hit else None
 # The fixed policy against every other filter: who has the lowest mean RMS.
 fixed={}
 for c in [x for x,_,_ in CELLS]:
  e={}
  for m in ('rms_z','rms_h'):
   e[m]={f:1000*mean(c,f,'heave_mm' if m=='rms_z' else 'horizontal_mm')/1000 for f in OPERATIONAL}
   e[m+'_best']=min(OPERATIONAL,key=lambda f:e[m][f])
   for f in ('CSM','TCSM'):e[m+'_'+f+'_vs_A41']=pct(c,f,'A41@C',m);e[m+'_'+f+'_vs_legacy']=pct(c,f,f+'/cap0',m)
   e[m+'_TCSM_vs_T41']=pct(c,'TCSM','T41@C',m);e[m+'_TCSM_vs_CSM']=pct(c,'TCSM','CSM',m)
  e['heading_deg']={f:mean(c,f) for f in OPERATIONAL+('CSM/cap0','TCSM/cap0')}
  e['heading_over30']={f:int(mean(c,f,'heading_over30')) for f in OPERATIONAL+('CSM/cap0','TCSM/cap0')}
  fixed[c]=e
 inter={}
 for c in INTERVENTION['cells']:
  g=iv[iv.cell==c];v=g.pivot(index='seed',columns='filter',values='rms_z');k=g[g['filter']=='CSM:noconf/ext-true/cap0/sf2']
  e=dict(known_coherence=float((k.cpl_coh_xz+k.cpl_coh_yz).mean()))
  for f,b in (('CSM/cap0/sf2','CSM/real/cap0/sf2'),('TCSM/cap0/sf2','TCSM/real/cap0/sf2')):
   n,d,lo,hi,_=ci(100*(v[f]/v[b]-1));e[f]=dict(heave_vs_real_pct=d,ci_low=lo,ci_high=hi)
  for f in INTERVENTION['filters']:
   x=g[g['filter']==f].yaw_rms_deg;e.setdefault(f,{}).update(heading_mean_deg=float(x.mean()),heading_over30=int((x>HEADING_FAILURE_DEG).sum()))
  inter[c]=e
 analysis=dict(analysis_sha256=sha(__file__),protocol_sha256=sha(out/'heave_sweep_protocol.json'),
  reproduction=reproduction(out,df),cap_binding=cap,legacy_heave_vs_real_only=harm,heading_onset=onset,
  control_heading_max_deg={f:max(mean(c,f) for c in DOSE) for f in ('A41@C','T41@C','CSM/real/cap0','TCSM/real/cap0')},
  intervention_sf2=inter,fixed_policy=fixed,
  spreading_at_015={f:{c:mean(c,f) for c in ('h0.15_broad','h0.15','h0.15_long')} for f in ('CSM/cap0','TCSM/cap0','CSM','TCSM')})
 (out/'heave_sweep_analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
 with pd.option_context('display.width',220):print(s.pivot(index='cell',columns='filter',values='heave_mm').round(2).to_string())
 print(json.dumps({k:analysis[k] for k in ('reproduction','heading_onset')},indent=1))

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['run','analyse'])
 ap.add_argument('--results',type=Path,default=OUT);ap.add_argument('--jobs',type=int,default=3)
 a=ap.parse_args();{'run':run,'analyse':analyse}[a.stage](a)
