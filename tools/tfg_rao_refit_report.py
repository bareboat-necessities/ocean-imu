#!/usr/bin/env python3
"""Analyse the TFG RAO core-coefficient refit: response surface, boundaries,
coefficient interactions, sea-state regimes, per-case consistency and the
predeclared promotion checks.

Reads only what the refit phases wrote.  The holdout sections appear only once
the sealed holdout has actually been executed.
"""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import tfg_rao_core_refit as R
from model_mismatch_ablation import RECORDS
HS={r.filename:r.hs_m for r in RECORDS}
KEYS,METRICS,W,PROTOCOL=R.KEYS,R.METRICS,R.W,R.PROTOCOL
PRIMARY=PROTOCOL['promotion_checks']['primary_metrics']
SECONDARY=PROTOCOL['promotion_checks']['secondary_metrics']
SHORT={'TFG_TAU_COEFF':'c_tau','TFG_SIGMA_COEFF':'c_sigma','TFG_R_S_X_FACTOR':'k_Sx','TFG_R_S_Y_FACTOR':'k_Sy'}

def composite(row,base):
 """Weighted log metric ratio for one paired case, or None if incomplete."""
 out=0.
 for m in METRICS:
  v,w=row['metrics'].get(m),base['metrics'].get(m)
  if v is None or w is None or not math.isfinite(v) or not math.isfinite(w) or w<=0 or v<0:return None
  out+=W[m]*math.log(max(v,1e-12)/w)
 return out/sum(W.values())

def per_case(runs):
 by=R.cases(runs);base=by['baseline'];out={}
 for name,cs in by.items():
  out[name]={k:composite(r,base[k]) for k,r in cs.items() if k in base}
 return by,base,out

def restricted(runs,keep):
 """Rank using only the cases a predicate keeps (regime and seed splits)."""
 sub=[r for r in runs if keep(r)]
 return R.rank(sub) if any(r['config']=='baseline' for r in sub) else []

def marginals(ranked,axis,fixed):
 """One-dimensional profile along one axis with the others held at `fixed`."""
 out=[]
 for r in ranked:
  c=r['coefficients']
  if all(abs(c[k]-v)<1e-9 for k,v in fixed.items()):out.append((c[axis],r))
 return sorted(out,key=lambda z:z[0])

def interactions(ranked):
 """Least-squares response surface in log-coefficients over the joint arms.

 Reports the pairwise interaction terms explicitly: the three coefficients act
 through tau, sigma_aw and r_S jointly, so a purely additive fit would hide the
 part of the surface the refinement has to search."""
 pts=[(r,[math.log(r['coefficients'][k]) for k in KEYS[:3]]) for r in ranked
      if r['eligible'] and math.isfinite(r['score'])]
 if len(pts)<12:return None
 y=np.array([r['raw_score'] for r,_ in pts])
 L=np.array([l for _,l in pts]);t,s,x=L[:,0],L[:,1],L[:,2]
 cols={'1':np.ones_like(t),'t':t,'s':s,'x':x,'t2':t*t,'s2':s*s,'x2':x*x,'ts':t*s,'tx':t*x,'sx':s*x}
 A=np.column_stack([cols[k] for k in cols])
 beta,*_=np.linalg.lstsq(A,y,rcond=None)
 resid=y-A@beta;ss=float(((y-y.mean())**2).sum())
 add=np.column_stack([cols[k] for k in ('1','t','s','x','t2','s2','x2')])
 badd,*_=np.linalg.lstsq(add,y,rcond=None);radd=y-add@badd
 return {'terms':dict(zip(cols,map(float,beta))),'n':len(pts),
         'r2':float(1-(resid**2).sum()/ss) if ss>0 else None,
         'r2_additive_only':float(1-(radd**2).sum()/ss) if ss>0 else None,
         'interaction_share_of_residual_variance':
            float(1-(resid**2).sum()/max((radd**2).sum(),1e-30))}

def consistency(comp,base_name='baseline'):
 """Per-case sign agreement of the composite improvement."""
 out={}
 for name,cs in comp.items():
  vals=[v for v in cs.values() if v is not None]
  if not vals:continue
  out[name]={'cases':len(vals),'improved_cases':sum(v<0 for v in vals),
             'unanimous':all(v<0 for v in vals) or all(v>0 for v in vals),
             'median':float(np.median(vals)),'worst_case':max(vals),'best_case':min(vals)}
 return out

def metric_table(runs,a,b):
 """Paired per-metric pooled ratio and per-record win counts for a vs b."""
 by=R.cases(runs);A,B=by[a],by[b];out={}
 for m in METRICS:
  logs=[];wins={};tot={}
  for k,row in A.items():
   o=B.get(k)
   if o is None:continue
   v,w=row['metrics'].get(m),o['metrics'].get(m)
   if v is None or w is None or not math.isfinite(v) or not math.isfinite(w) or w<=0:continue
   logs.append(math.log(max(v,1e-12)/w))
   rec=k[0];tot[rec]=tot.get(rec,0)+1;wins[rec]=wins.get(rec,0)+(v<w)
  if not logs:out[m]=None;continue
  out[m]={'pooled_pct_change':float((math.exp(float(np.mean(logs)))-1)*100),
          'cases':len(logs),'cases_improved':int(sum(l<0 for l in logs)),
          'records_majority_improved':int(sum(wins.get(r,0)*2>tot[r] for r in tot)),
          'records':len(tot),'per_record_pct':{r:float((math.exp(float(np.mean(
             [l for l,(kk) in zip(logs,[k for k in A if k in B])if kk[0]==r])))-1)*100) for r in tot}}
 return out

def fmt(v,n=3):
 return 'n/a' if v is None else ('inf' if isinstance(v,float) and not math.isfinite(v) else f'{v:.{n}f}')

def surface_md(title,ranked,limit=None):
 lines=[f'### {title}','',
        '| rank | arm | c_tau | c_sigma | k_Sx | k_Sy | score | '+' | '.join(f'{m} x' for m in METRICS)+' | gate fails | viol | eligible |',
        '|---:|---|---:|---:|---:|---:|---:|'+'---:|'*len(METRICS)+'---:|---:|---|']
 for i,r in enumerate(ranked if limit is None else ranked[:limit],1):
  c=r['coefficients']
  lines.append('| '+' | '.join([str(i),r['config']]+[fmt(c[k],3) for k in KEYS]+[fmt(r['score'],5)]+
               [fmt(r['ratios'][m],4) for m in METRICS]+
               [str(r['gate_failures']),str(r['violations']),'yes' if r['eligible'] else 'NO'])+' |')
 return lines+['']

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_rao_core_refit')
 a=p.parse_args();o=a.output_dir.resolve()
 doc=['# TFG RAO core-coefficient refit','',
      'Pinned 28-ft vessel-RAO records '+str(PROTOCOL['record_indices'])+'; '
      f"Q_bg fixed at {PROTOCOL['fixed_gyro_bias_rw_var']}.  Ratios are the candidate over the",
      'paired baseline case, so values below 1 are improvements.','']
 J={'protocol':PROTOCOL}
 phases=[t for t in ('screen','extend','refine','holdout') if (o/t/'runs.json').exists()]
 doc+=['Executed phases: '+', '.join(phases) if phases else 'No phase has produced runs yet.','']
 train=[t for t in ('screen','extend') if t in phases]

 if train:
  runs=R.load(o,train);ranked=R.rank(runs)
  (o/'screen-ranking.json').write_text(json.dumps(ranked,indent=2)+'\n')
  by,base,comp=per_case(runs)
  doc+=['## 1. Training response surface','',
        f'{len(ranked)} arms x {len(base)} paired cases '
        f"(8 records x {len(PROTOCOL['training_seeds'])} training seeds).",'']
  doc+=surface_md('Full ranking (all arms retained, including ineligible ones)',ranked)
  best=next(r for r in ranked if r['eligible'])
  bc=best['coefficients'];axes=R.swept()
  doc+=['## 2. Boundary contacts','',
        '| axis | best | swept min | swept max | on boundary |','|---|---:|---:|---:|---|']
  hits=[]
  for k in KEYS[:3]:
   lo,hi=axes[k][0],axes[k][-1];hit='low' if abs(bc[k]-lo)<1e-9 else ('high' if abs(bc[k]-hi)<1e-9 else '')
   if hit:hits.append({'axis':k,'side':hit,'value':bc[k]})
   doc.append(f'| {SHORT[k]} | {bc[k]:.4g} | {lo:.4g} | {hi:.4g} | {hit or "no"} |')
  doc+=['',('Boundary contact present: the screen optimum sits on a search face, so that axis '
            'is extended before refinement.') if hits else
           'No axis optimum sits on a search boundary; the optimum is interior to the screen.','']
  J['boundary_contacts']=hits
  doc+=['## 3. One-dimensional marginals through the pre-RAO point','']
  for k in KEYS[:3]:
   fixed={j:v for j,v in zip(KEYS,R.LEGACY) if j!=k and not (k==KEYS[2] and j==KEYS[3])}
   if k==KEYS[2]:fixed={KEYS[0]:R.LEGACY[0],KEYS[1]:R.LEGACY[1]}
   prof=marginals(ranked,k,fixed)
   if len(prof)<3:continue
   doc+=[f'{SHORT[k]}: '+'  '.join(f'{v:.4g}->{r["raw_score"]:+.4f}' for v,r in prof),'']
  inter=interactions(ranked);J['interactions']=inter
  if inter:
   doc+=['## 4. Coefficient interactions','',
         'Quadratic least squares of the raw score on log c_tau, log c_sigma, log k_S.','',
         '| term | coefficient |','|---|---:|']
   doc+=[f'| {t} | {v:+.5f} |' for t,v in inter['terms'].items()]
   doc+=['',f"Fit R2 {fmt(inter['r2'],4)}; additive-only R2 {fmt(inter['r2_additive_only'],4)}; "
          f"the cross terms remove {inter['interaction_share_of_residual_variance']*100:.1f}% of the "
          'additive model residual variance.',
          '','Non-zero cross terms mean the axes cannot be optimised one at a time; the refinement '
          'is therefore a joint grid, not three independent sweeps.','']
  doc+=['## 5. Sea-state regimes','']
  lows=[f for f,h in HS.items() if h<=1.5];highs=[f for f,h in HS.items() if h>1.5]
  reg={}
  for label,keep in (('low Hs (0.27, 1.50 m)',lambda r:r['input'] in lows),
                     ('high Hs (4.00, 8.50 m)',lambda r:r['input'] in highs)):
   rr=[x for x in restricted(runs,keep) if x['eligible']]
   if not rr:continue
   reg[label]={'best':rr[0]['config'],'coefficients':rr[0]['coefficients'],'score':rr[0]['score']}
   doc.append(f"- {label}: best arm `{rr[0]['config']}` "+
              ', '.join(f'{SHORT[k]}={rr[0]["coefficients"][k]:.4g}' for k in KEYS[:3])+
              f" (score {rr[0]['score']:+.4f})")
  J['regimes']=reg
  if len(reg)==2:
   a1,b1=(v['coefficients'] for v in reg.values())
   spread={SHORT[k]:(b1[k]/a1[k] if a1[k] else None) for k in KEYS[:3]}
   J['regime_ratio']=spread
   doc+=['',f'High/low optimum ratio per axis: '+', '.join(f'{k}={fmt(v,3)}' for k,v in spread.items()),
         '',('The two regimes prefer materially different coefficients, which is a question for the '
             'adaptation law rather than for a single global constant.')
          if any(v and (v>1.25 or v<0.8) for v in spread.values()) else
          ('The two regimes prefer the same basin, so one global constant per axis remains the right '
           'shape for these records.'),'']
  cons=consistency(comp);J['training_consistency']={k:cons[k] for k in list(cons)[:400]}
  doc+=['## 6. Per-case consistency of the leading training arms','',
        '| arm | cases | improved | unanimous | median log ratio | worst case |','|---|---:|---:|---|---:|---:|']
  for r in ranked[:12]:
   c=cons.get(r['config'])
   if not c:continue
   doc.append(f"| {r['config']} | {c['cases']} | {c['improved_cases']} | "
              f"{'yes' if c['unanimous'] else 'no'} | {c['median']:+.5f} | {c['worst_case']:+.5f} |")
  doc+=['']
  J['training_best']={'config':best['config'],'coefficients':bc,'score':best['score'],
                      'ratios':best['ratios']}

 if 'refine' in phases:
  rruns=R.load(o,['refine']);rranked=R.rank(rruns)
  (o/'refine-ranking.json').write_text(json.dumps(rranked,indent=2)+'\n')
  _,rbase,rcomp=per_case(rruns)
  doc+=['## 7. Refinement (independent refinement seeds '+
        ', '.join(PROTOCOL['refinement_seeds'])+')','']
  doc+=surface_md('Refinement ranking',rranked)
  iso=[r for r in rranked if r['eligible'] and
       abs(r['coefficients'][KEYS[2]]-r['coefficients'][KEYS[3]])<1e-9]
  ani=[r for r in rranked if r['eligible'] and
       abs(r['coefficients'][KEYS[2]]-r['coefficients'][KEYS[3]])>=1e-9]
  doc+=['### X/Y anisotropy','']
  if iso and ani:
   doc+=[f"Best isotropic arm `{iso[0]['config']}` scores {iso[0]['score']:+.5f}; "
         f"best anisotropic arm `{ani[0]['config']}` scores {ani[0]['score']:+.5f} "
         f"(k_Sx={ani[0]['coefficients'][KEYS[2]]:.4g}, k_Sy={ani[0]['coefficients'][KEYS[3]]:.4g}).",'',
         ('Anisotropy wins on the refinement seeds.' if ani[0]['score']<iso[0]['score'] else
          'Anisotropy does not beat the isotropic arm, so k_Sx = k_Sy is kept.'),'']
   J['anisotropy']={'best_isotropic':iso[0]['config'],'isotropic_score':iso[0]['score'],
                    'best_anisotropic':ani[0]['config'],'anisotropic_score':ani[0]['score'],
                    'anisotropy_wins':bool(ani[0]['score']<iso[0]['score'])}
  rc=consistency(rcomp)
  doc+=['| arm | cases | improved | unanimous | median log ratio | worst case |','|---|---:|---:|---|---:|---:|']
  for r in rranked[:12]:
   c=rc.get(r['config'])
   if c:doc.append(f"| {r['config']} | {c['cases']} | {c['improved_cases']} | "
                   f"{'yes' if c['unanimous'] else 'no'} | {c['median']:+.5f} | {c['worst_case']:+.5f} |")
  doc+=['']

 if (o/'candidate.json').exists():
  cand=json.loads((o/'candidate.json').read_text());J['candidate']=cand
  doc+=['## 8. Frozen candidate','',
        '| coefficient | pre-RAO baseline | candidate |','|---|---:|---:|']
  for k,l in zip(KEYS,R.LEGACY):
   doc.append(f'| {SHORT[k]} | {l:.4g} | {cand["coefficients"][k]:.4g} |')
  doc+=['',f"Frozen from `{cand['config']}` on training and refinement seeds only; "
        f"coefficient SHA-256 `{cand['coefficient_sha256'][:16]}...`.",'']

 if 'holdout' in phases:
  hruns=R.load(o,['holdout']);hranked=R.rank(hruns)
  doc+=['## 9. Sealed holdout (seeds '+', '.join(PROTOCOL['sealed_holdout_seeds'])+')','']
  doc+=surface_md('Holdout ranking',hranked)
  tab=metric_table(hruns,'candidate','baseline');J['holdout_vs_baseline']=tab
  doc+=['### Candidate against the current TFG baseline','',
        '| metric | pooled % change | cases improved | records with a majority improved |',
        '|---|---:|---:|---:|']
  for m in METRICS:
   t=tab.get(m)
   doc.append(f"| {m} | {t['pooled_pct_change']:+.3f} | {t['cases_improved']}/{t['cases']} | "
              f"{t['records_majority_improved']}/{t['records']} |" if t else f'| {m} | n/a | n/a | n/a |')
  doc+=['','Per-record % change (negative is better):','',
        '| metric | '+' | '.join(f'Hs {HS[f]:g} {f.split("_")[2]}' for f in sorted(HS,key=lambda z:(HS[z],z)))+' |',
        '|---|'+'---:|'*len(HS)]
  for m in METRICS:
   t=tab.get(m)
   if not t:continue
   doc.append(f'| {m} | '+' | '.join(f"{t['per_record_pct'].get(f,float('nan')):+.2f}"
              for f in sorted(HS,key=lambda z:(HS[z],z)))+' |')
  doc+=['']
  ou=o/'holdout-ou3'/'runs.json'
  if ou.exists():
   both=hruns+json.loads(ou.read_text())
   ot=metric_table(both,'candidate','ou3_shipping');J['holdout_vs_ou3']=ot
   bt=metric_table(both,'baseline','ou3_shipping');J['baseline_vs_ou3']=bt
   doc+=['### Against shipping OU-III on identical holdout histories','',
         '| metric | candidate vs OU-III % | baseline vs OU-III % |','|---|---:|---:|']
   for m in METRICS:
    doc.append(f"| {m} | {fmt(ot[m]['pooled_pct_change'],3) if ot.get(m) else 'n/a'} | "
               f"{fmt(bt[m]['pooled_pct_change'],3) if bt.get(m) else 'n/a'} |")
   doc+=['']
  pc=PROTOCOL['promotion_checks'];by=R.cases(hruns);checks=[]
  cand_rows,base_rows=by.get('candidate',{}),by.get('baseline',{})
  complete=all(all(m in r['metrics'] and math.isfinite(r['metrics'][m]) for m in METRICS)
               for r in list(cand_rows.values())+list(base_rows.values())) and \
           len(cand_rows)==len(base_rows)==8*len(PROTOCOL['sealed_holdout_seeds'])
  checks.append(('all cases complete and finite',complete,f'{len(cand_rows)} candidate cases'))
  newfail=sum(1 for k,r in cand_rows.items() if r['exit_code']!=0 and base_rows[k]['exit_code']==0)
  checks.append(('no new gate-failure cases',newfail<=pc['new_gate_failure_cases_vs_tfg_baseline'],
                 f'{newfail} new failures'))
  extra=sum(len(r['violations']) for r in cand_rows.values())-sum(len(r['violations']) for r in base_rows.values())
  checks.append(('no additional gate violations',extra<=pc['additional_gate_violations_vs_tfg_baseline'],
                 f'{extra:+d} violations'))
  for m in PRIMARY:
   t=tab.get(m)
   checks.append((f'{m} pooled regression <= {pc["maximum_pooled_primary_metric_regression_pct"]}%',
                  bool(t and t['pooled_pct_change']<=pc['maximum_pooled_primary_metric_regression_pct']),
                  fmt(t and t['pooled_pct_change'],3)+'%'))
   checks.append((f'{m} improved on >= {pc["minimum_record_wins_for_each_primary_metric"]} records',
                  bool(t and t['records_majority_improved']>=pc['minimum_record_wins_for_each_primary_metric']),
                  f"{t['records_majority_improved']}/8" if t else 'n/a'))
  for m in SECONDARY:
   t=tab.get(m)
   checks.append((f'{m} pooled regression <= {pc["maximum_pooled_secondary_metric_regression_pct"]}%',
                  bool(t and t['pooled_pct_change']<=pc['maximum_pooled_secondary_metric_regression_pct']),
                  fmt(t and t['pooled_pct_change'],3)+'%'))
  cscore=next((r for r in hranked if r['config']=='candidate'),None)
  comp_pct=(math.exp(cscore['raw_score'])-1)*100 if cscore and math.isfinite(cscore['raw_score']) else None
  checks.append((f'balanced composite improves by >= {pc["minimum_balanced_composite_improvement_pct"]}%',
                 bool(comp_pct is not None and -comp_pct>=pc['minimum_balanced_composite_improvement_pct']),
                 fmt(comp_pct,3)+'%'))
  doc+=['### Predeclared promotion checks','','| check | result | observed |','|---|---|---:|']
  doc+=[f'| {n} | {"PASS" if ok else "FAIL"} | {obs} |' for n,ok,obs in checks]
  verdict=all(ok for _,ok,_ in checks)
  J['promotion']={'checks':[{'check':n,'pass':bool(ok),'observed':obs} for n,ok,obs in checks],
                  'promote':bool(verdict)}
  doc+=['',f'**Verdict: {"PROMOTE" if verdict else "DO NOT PROMOTE"}** '
        '— the sealed holdout decides; a candidate that only won training is not promoted.','']

 (o/'report.md').write_text('\n'.join(doc)+'\n')
 (o/'analysis.json').write_text(json.dumps(J,indent=2,default=str)+'\n')
 print('\n'.join(doc[:60]));print(f'\n[wrote {o/"report.md"} and {o/"analysis.json"}]')
if __name__=='__main__':main()
