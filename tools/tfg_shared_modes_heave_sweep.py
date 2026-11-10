#!/usr/bin/env python3
"""Exploratory reduced-heave dose-response for the shared-mode filters.

The primary study finds shared-mode heading collapse at heave RAO x0.15 and
none at x1. This sweep locates the onset on the primary seeds and tests two
predictions of the working hypothesis that the vertical-horizontal quadrature
prior, not the TFG geometry, drives the collapse:

  Dose       shared heading error grows as the heave RAO shrinks, while real-only
             and independent-axis controls stay near the nominal level;
  Spreading  at x0.15, broad spreading (lower vertical-horizontal coherence) weakens
             the collapse and long-crested seas strengthen it.

Two later predictions were registered from the first analysis (EXTENSION and
INTERVENTION below): the positive-semidefinite cap of the excitation mapping
binds below a computable heave factor, and raising only the horizontal prior
moves that boundary.

Known coupling is a diagnostic input only. No estimator is retuned and the
results are outside the primary multiplicity family. The x1 and default x0.15
cells repeat registered seeds, so the overlapping per-run scores must
reproduce the published tables exactly.

Run:        python3 tools/tfg_shared_modes_heave_sweep.py run --jobs 4
Reanalyse:  python3 tools/tfg_shared_modes_heave_sweep.py analyse
"""
from __future__ import annotations
import argparse,csv,io,json,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from tfg_shared_modes_study import OUT,BINARY,BASE,clean_env,sha,ci

FILTERS=('A41@C','CSM','TCSM','CSM/real','TCSM/real','CSM:noconf/ext-true','TCSM:noconf/ext-true',
         'TCSM:noconf/ext-true/real','TCSM/noreset')
SEEDS=list(range(1001,1013))
# (cell, heave RAO factor, spreading exponent s; None keeps the default s=4)
CELLS=(('h1.00',1.,None),('h0.80',.8,None),('h0.70',.7,None),('h0.60',.6,None),('h0.50',.5,None),('h0.40',.4,None),('h0.30',.3,None),('h0.25',.25,None),
       ('h0.20',.2,None),('h0.15',.15,None),('h0.15_broad',.15,1),('h0.15_long',.15,10000))
DOSE=[c for c,_,s in CELLS if s is None]
# Registered after the first analysis: the known prior's applied multiple
# coherence is 0.334 at x1 and pinned at the 0.999 PSD cap for every factor
# <= 0.5. With the cross term scaled by the vertical prior variance, the
# uncapped coherence grows like 0.334 / k^2, so the cap binds for
# k < sqrt(0.334 / 0.999) = 0.58. These cells test whether shared-prior heave
# harm begins only where the cap binds.
EXTENSION=('h0.80','h0.70','h0.60')
CAP=.999
# Registered after the extension confirmed the 1/k^2 law (0.523, 0.685, 0.936
# at k = 0.8, 0.7, 0.6): an intervention on the cap alone. The wrapper's
# horizontal prior is S_factor times the vertical one (default 1); S_factor 2
# divides the applied coherence by four, so the cap should bind only for
# k < sqrt(0.334 / (4 * 0.999)) = 0.29. Same seeds, cells across the boundary.
INTERVENTION=dict(s_factor=2.,cells=('h1.00','h0.50','h0.30','h0.25','h0.20'),
 filters=('A41@C/sf2','CSM/sf2','TCSM/sf2','CSM/real/sf2','TCSM/real/sf2','CSM:noconf/ext-true/sf2','TCSM:noconf/ext-true/sf2'))
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
  role='exploratory dose-response diagnosis; no estimator retuning; outside the primary multiplicity family',
  extension=dict(cells=EXTENSION,registered='after the first analysis of the other cells',
   prediction='cap: the 0.999 PSD cap of the known prior binds for k < 0.58 (from coherence 0.334 at x1); shared-prior heave harm relative to real-only coupling begins only where the cap binds'),
  intervention=dict(INTERVENTION,registered='after the extension analysis',
   prediction='intervention: with S_factor 2 the known prior reaches the cap only for k < 0.29; shared heave harm relative to real-only coupling and the heading collapse are absent at k >= 0.3 and return at k <= 0.25'),
  predictions={'dose':'shared heading error increases as heave RAO decreases; real-only and independent controls stay nominal',
               'spreading':'at x0.15 broad spreading weakens and long-crested seas strengthen the shared heading collapse'})

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
 """Summaries, paired contrasts against controls, and the onset of collapse."""
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
 paired=[]
 for cell,g in df.groupby('cell'):
  for metric in ('yaw_rms_deg','rms_z','rms_h'):
   v=g.pivot(index='seed',columns='filter',values=metric)
   for a,b in (('CSM','A41@C'),('TCSM','A41@C'),('TCSM','CSM'),('CSM','CSM/real'),('TCSM','TCSM/real'),
               ('TCSM/noreset','TCSM'),('TCSM:noconf/ext-true','CSM:noconf/ext-true'),
               ('TCSM:noconf/ext-true','TCSM:noconf/ext-true/real')):
    n,d,lo,hi,p=ci(v[a]-v[b]);paired.append(dict(cell=cell,filter=a,baseline=b,metric=metric,n=n,delta=d,ci_low=lo,ci_high=hi,p=p))
 # Spreading contrasts at x0.15 pair each seed across the three sea cells.
 for f in FILTERS:
  for other in ('h0.15_broad','h0.15_long'):
   v=df[df['filter']==f].pivot(index='seed',columns='cell',values='yaw_rms_deg')
   n,d,lo,hi,p=ci(v[other]-v['h0.15']);paired.append(dict(cell=other+'-h0.15',filter=f,baseline=f,metric='yaw_rms_deg',n=n,delta=d,ci_low=lo,ci_high=hi,p=p))
 pd.DataFrame(paired).to_csv(out/'heave_sweep_paired.csv',index=False,float_format='%.9g')
 s=summary.reset_index()
 # Onset: first heave factor (default spreading) whose paired shared-minus-
 # independent heading interval excludes zero by more than one degree.
 pr=pd.DataFrame(paired);onset={}
 for f in ('CSM','TCSM'):
  hit=[c for c in DOSE if pr[(pr.cell==c)&(pr['filter']==f)&(pr.baseline=='A41@C')&(pr.metric=='yaw_rms_deg')].ci_low.iloc[0]>1.]
  onset[f]=hit[0] if hit else None
 def mean(cell,f):return float(s[(s.cell==cell)&(s['filter']==f)].heading_deg.iloc[0])
 shared_mono={f:bool(np.all(np.diff([mean(c,f) for c in DOSE])>0)) for f in ('CSM','TCSM')}
 controls_max={f:max(mean(c,f) for c in DOSE) for f in ('A41@C','CSM/real','TCSM/real')}
 mc=df.assign(mc=df.cpl_coh_xz+df.cpl_coh_yz).groupby(['cell','filter']).mc.mean()
 cap={c:dict(known=float(mc[(c,'CSM:noconf/ext-true')]),estimated_qmekf=float(mc[(c,'CSM')]),estimated_tfg=float(mc[(c,'TCSM')]),
            known_at_cap=bool(mc[(c,'CSM:noconf/ext-true')]>=CAP-5e-4)) for c in DOSE}
 harm={}
 for f,b in (('CSM','CSM/real'),('TCSM','TCSM/real')):
  hit=[]
  for c in DOSE:
   v=df[df.cell==c].pivot(index='seed',columns='filter',values='rms_z')
   n,d,lo,hi,_=ci(100*(v[f]/v[b]-1));harm.setdefault(f,{})[c]=dict(percent=d,ci_low=lo,ci_high=hi)
   if lo>0:hit.append(c)
  harm[f]['first_harm_cell']=hit[0] if hit else None
 inter={}
 for c in INTERVENTION['cells']:
  g=iv[iv.cell==c];v=g.pivot(index='seed',columns='filter',values='rms_z');k=g[g['filter']=='CSM:noconf/ext-true/sf2']
  e=dict(known_coherence=float((k.cpl_coh_xz+k.cpl_coh_yz).mean()))
  for f,b in (('CSM/sf2','CSM/real/sf2'),('TCSM/sf2','TCSM/real/sf2')):
   n,d,lo,hi,_=ci(100*(v[f]/v[b]-1));e[f]=dict(heave_vs_real_pct=d,ci_low=lo,ci_high=hi)
  for f in ('A41@C/sf2','CSM/sf2','TCSM/sf2','CSM/real/sf2','TCSM/real/sf2','CSM:noconf/ext-true/sf2','TCSM:noconf/ext-true/sf2'):
   x=g[g['filter']==f].yaw_rms_deg;e.setdefault(f,{}).update(heading_mean_deg=float(x.mean()),heading_over30=int((x>HEADING_FAILURE_DEG).sum()))
  inter[c]=e
 analysis=dict(analysis_sha256=sha(__file__),cap_binding=cap,heave_vs_real_only=harm,intervention_sf2=inter,protocol_sha256=sha(out/'heave_sweep_protocol.json'),
  reproduction=reproduction(out,df),onset_cell=onset,shared_heading_strictly_increasing=shared_mono,
  control_heading_max_deg=controls_max,
  spreading_at_015={f:{c:mean(c,f) for c in ('h0.15_broad','h0.15','h0.15_long')} for f in ('CSM','TCSM')})
 (out/'heave_sweep_analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
 with pd.option_context('display.width',200):print(summary.to_string())
 print(json.dumps(analysis,indent=1))

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['run','analyse'])
 ap.add_argument('--results',type=Path,default=OUT);ap.add_argument('--jobs',type=int,default=4)
 a=ap.parse_args();{'run':run,'analyse':analyse}[a.stage](a)
