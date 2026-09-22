#!/usr/bin/env python3
"""Focused TFG RAO core-coefficient refit; pre-RAO coefficients are baselines only.

Each phase is an independent operation using the seed set predeclared in
tools/tfg_rao_refit_protocol.json.  Screening, boundary extension, refinement
and candidate freezing never read a holdout seed, and the holdout phase never
selects.  Failed and non-finite arms are retained as evidence and marked
ineligible rather than dropped.
"""
from __future__ import annotations
import argparse,hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'tools/rao_parameter_tuning.py'
PROTOCOL=json.loads((ROOT/'tools/tfg_rao_refit_protocol.json').read_text())
SEEDS={'screen':PROTOCOL['training_seeds'],'extend':PROTOCOL['training_seeds'],
       'refine':PROTOCOL['refinement_seeds'],'holdout':PROTOCOL['sealed_holdout_seeds']}
RECORDS=','.join(str(i) for i in PROTOCOL['record_indices'])
BASEENV={'TFG_TUNING':'adaptive','TFG_AW_COV_SYNC':'1',
         'SF_GYRO_BIAS_RW_VAR':repr(PROTOCOL['fixed_gyro_bias_rw_var'])}
W=PROTOCOL['physical_metric_weights']; METRICS=tuple(W)
KEYS=tuple(PROTOCOL['coefficient_keys'])
# Pre-RAO shipping coefficients.  Baseline only: the sweeps that produced them
# predate the pinned 28-ft vessel-RAO dataset.
LEGACY=(1.0,0.8,1.15,1.15)

def arm(name,tau=LEGACY[0],sigma=LEGACY[1],x=LEGACY[2],y=LEGACY[3]):
 return {'name':name,'env':dict(BASEENV,**dict(zip(KEYS,(str(tau),str(sigma),str(x),str(y)))))}
def coeffs(env):
 return tuple(float(env.get(k,d)) for k,d in zip(KEYS,LEGACY))
def dedupe(arms):
 seen=set();out=[]
 for a in arms:
  k=coeffs(a['env'])
  if k in seen:continue
  seen.add(k);out.append(a)
 return out
def worst(row):
 return max([v['ratio'] for v in row['violations'] if v['ratio'] is not None],default=1.)

def cases(runs):
 """Per-configuration paired cases keyed by (record, seed)."""
 by={}
 for r in runs:by.setdefault(r['config'],{})[(r['input'],r['seed'])]=r
 return by

def rank(runs):
 """Predeclared score: the equal-case mean of weighted log metric ratios against
 the paired baseline case, plus the gate penalty.  Every record and every seed
 carries equal weight.  An arm with a missing or non-finite metric is retained
 and marked ineligible; a missing metric is never read as zero."""
 by=cases(runs);base=by['baseline'];ans=[]
 for name,cs in by.items():
  logs=[];ratios={m:[] for m in METRICS};ok=len(cs)==len(base)
  for k,row in cs.items():
   b=base.get(k)
   if b is None:ok=False;continue
   vals=[(row['metrics'].get(m),b['metrics'].get(m)) for m in METRICS]
   if any(v is None or w is None or not math.isfinite(v) or not math.isfinite(w) or w<=0 or v<0
          for v,w in vals):ok=False;continue
   for m,(v,w) in zip(METRICS,vals):ratios[m].append(max(v,1e-12)/w)
   logs.append(sum(W[m]*math.log(max(v,1e-12)/w) for m,(v,w) in zip(METRICS,vals))/sum(W.values()))
  viol=sum(len(r['violations']) for r in cs.values())
  bviol=sum(len(r['violations']) for r in base.values())
  penalty=.025*max(0,viol-bviol)+.025*max(0.,max(map(worst,cs.values()))-max(1.,max(map(worst,base.values()))))
  raw=sum(logs)/len(logs) if logs else math.inf
  good=bool(ok and logs)
  ans.append(dict(config=name,coefficients=dict(zip(KEYS,coeffs(next(iter(cs.values()))['env']))),
                  eligible=good,cases=len(logs),expected_cases=len(base),
                  score=raw+penalty if good else math.inf,raw_score=raw,penalty=penalty,
                  gate_failures=sum(r['exit_code']!=0 for r in cs.values()),violations=viol,
                  ratios={m:(math.exp(sum(map(math.log,v))/len(v)) if v else None) for m,v in ratios.items()}))
 return sorted(ans,key=lambda z:(z['score'],z['config']))

def run(tag,cfgs,seeds,out,jobs,family='TFG'):
 d=out/tag;d.mkdir(parents=True,exist_ok=True);p=d/'configs.json'
 p.write_text(json.dumps(cfgs,indent=2)+'\n')
 subprocess.run([sys.executable,str(RUN),'--family',family,'--configs',str(p),'--records',RECORDS,
                 '--seeds',','.join(seeds),'--jobs',str(jobs),'--output-dir',str(d)],cwd=ROOT,check=True)
 return json.loads((d/'runs.json').read_text())

def screen():
 """Deliberately far wider than the pre-RAO operating point."""
 out=[arm('baseline')]
 for v in (.45,.6,.75,.9,1.1,1.3,1.55,1.85):out.append(arm(f'tau_{v}',tau=v))
 for v in (.35,.5,.65,.95,1.15,1.4):out.append(arm(f'sigma_{v}',sigma=v))
 for v in (.45,.65,.85,1.0,1.3,1.55,1.85,2.2):out.append(arm(f'rsxy_{v}',x=v,y=v))
 for t in (.65,.9,1.15,1.45):
  for s in (.5,.75,1.0,1.25):
   for xy in (.65,1.0,1.35,1.75):out.append(arm(f'joint_t{t}_s{s}_xy{xy}',t,s,xy,xy))
 return dedupe(out)

def swept():
 """Per-axis values the screen actually visits, for boundary detection."""
 axes={k:set() for k in KEYS}
 for a in screen():
  for k,v in zip(KEYS,coeffs(a['env'])):axes[k].add(round(v,6))
 return {k:sorted(v) for k,v in axes.items()}

def extend(best):
 """Extend only the axes whose screen optimum sits on a search boundary."""
 t,s,x,_=best; axes=swept(); arms=[arm('baseline')]; hits=[]
 for k,v in zip(KEYS[:3],(t,s,x)):
  for side,edge in (('low',axes[k][0]),('high',axes[k][-1])):
   if abs(v-edge)>1e-9:continue
   hits.append({'axis':k,'side':side,'value':v})
   for f in ((.78,.6,.46) if side=='low' else (1.28,1.65,2.1)):
    c=dict(zip(KEYS[:3],(t,s,x)));c[k]=v*f
    arms.append(arm(f"ext_{k.split('_',1)[1].lower()}_{v*f:.4g}",
                    c[KEYS[0]],c[KEYS[1]],c[KEYS[2]],c[KEYS[2]]))
 return dedupe(arms),hits

def refine(best):
 """Finer joint grid inside the winning basin, then explicit X/Y anisotropy."""
 t,s,x,y=best; arms=[arm('baseline')]
 for tt in (t*.85,t,t*1.18):
  for ss in (s*.85,s,s*1.18):
   for xx in (x*.85,x,x*1.18):
    arms.append(arm(f'r_t{tt:.4g}_s{ss:.4g}_xy{xx:.4g}',tt,ss,xx,xx))
 for fx in (.8,.9,1.,1.12,1.25):
  for fy in (.8,.9,1.,1.12,1.25):
   arms.append(arm(f'aniso_x{x*fx:.4g}_y{y*fy:.4g}',t,s,x*fx,y*fy))
 return dedupe(arms)

def load(o,tags):
 rows=[];seen=set()
 for tag in tags:
  f=o/tag/'runs.json'
  if not f.exists():continue
  for r in json.loads(f.read_text()):
   key=(r['config'],r['input'],r['seed'])
   if key in seen:continue
   seen.add(key);rows.append(r)
 if not rows:raise SystemExit(f'no completed phase among {tags}')
 return rows

def best_of(o,tags):
 top=[r for r in rank(load(o,tags)) if r['eligible']]
 if not top:raise SystemExit('no eligible arm')
 return top[0]

def status(o,phase,**kw):
 (o/f'phase-status-{phase}.json').write_text(json.dumps(dict(phase=phase,**kw),indent=2)+'\n')

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_rao_core_refit')
 p.add_argument('--jobs',type=int,default=4)
 p.add_argument('--phase',choices=('screen','extend','refine','freeze','holdout'),default='screen')
 p.add_argument('--shard-index',type=int,default=0)
 p.add_argument('--shard-count',type=int,default=1)
 a=p.parse_args()
 if not 0<=a.shard_index<a.shard_count:p.error('invalid shard index/count')
 o=a.output_dir.resolve();o.mkdir(parents=True,exist_ok=True)

 if a.phase=='screen':
  cs=screen()[a.shard_index::a.shard_count]
  rows=run('screen',cs,SEEDS['screen'],o,a.jobs)
  if a.shard_count==1:(o/'screen-ranking.json').write_text(json.dumps(rank(rows),indent=2)+'\n')
  status(o,'screen',shard_index=a.shard_index,shard_count=a.shard_count,configurations=len(cs),
         records=len(PROTOCOL['record_indices']),seeds=SEEDS['screen'],holdout_evaluated=False)
 elif a.phase=='extend':
  top=best_of(o,['screen'])
  cs,hits=extend(tuple(top['coefficients'][k] for k in KEYS))
  (o/'boundary-contacts.json').write_text(json.dumps(
      {'screen_best':top['config'],'coefficients':top['coefficients'],'contacts':hits},indent=2)+'\n')
  if len(cs)<2:
   status(o,'extend',configurations=0,boundary_contacts=hits,seeds=SEEDS['extend'],
          holdout_evaluated=False);return
  rows=run('extend',cs,SEEDS['extend'],o,a.jobs)
  (o/'extend-ranking.json').write_text(json.dumps(rank(rows),indent=2)+'\n')
  status(o,'extend',configurations=len(cs),boundary_contacts=hits,seeds=SEEDS['extend'],
         holdout_evaluated=False)
 elif a.phase=='refine':
  top=best_of(o,['screen','extend'])
  cs=refine(tuple(top['coefficients'][k] for k in KEYS))
  (o/'refine-center.json').write_text(json.dumps(
      {'center_config':top['config'],'coefficients':top['coefficients']},indent=2)+'\n')
  rows=run('refine',cs,SEEDS['refine'],o,a.jobs)
  (o/'refine-ranking.json').write_text(json.dumps(rank(rows),indent=2)+'\n')
  status(o,'refine',configurations=len(cs),seeds=SEEDS['refine'],holdout_evaluated=False)
 elif a.phase=='freeze':
  top=[r for r in json.loads((o/'refine-ranking.json').read_text()) if r['eligible']][0]
  c=top['coefficients']
  (o/'candidate.json').write_text(json.dumps(
      {'config':top['config'],'coefficients':c,
       'coefficient_sha256':hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest(),
       'selected_on':'training and refinement seeds only','refinement_score':top['score'],
       'holdout_seeds_used_for_selection':False},indent=2)+'\n')
  status(o,'freeze',candidate=top['config'],coefficients=c,holdout_evaluated=False)
 else:
  cand=json.loads((o/'candidate.json').read_text());c=cand['coefficients']
  if hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()!=cand['coefficient_sha256']:
   raise SystemExit('frozen candidate hash mismatch')
  cs=dedupe([arm('baseline'),arm('candidate',*[c[k] for k in KEYS])])
  rows=run('holdout',cs,SEEDS['holdout'],o,a.jobs)
  (o/'holdout-ranking.json').write_text(json.dumps(rank(rows),indent=2)+'\n')
  run('holdout-ou3',[{'name':'ou3_shipping','env':{}}],SEEDS['holdout'],o,a.jobs,family='OU-III')
  status(o,'holdout',candidate=cand['config'],coefficients=c,seeds=SEEDS['holdout'],
         configurations=len(cs),comparators=PROTOCOL['holdout_comparators'],holdout_evaluated=True)
if __name__=='__main__':main()
