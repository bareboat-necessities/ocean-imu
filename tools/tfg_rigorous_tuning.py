#!/usr/bin/env python3
"""TFG-specific staged tuning with sealed holdout evaluation.

The campaign optimizes physical RMS endpoints, not tangent-coordinate covariance.
Candidate selection uses training/refinement draws only. The four holdout draws
are evaluated exactly once after the final candidate is frozen. Quality gates
are recorded but never changed or used to discard a replay.
"""
from __future__ import annotations
import argparse,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/'tools/rao_parameter_tuning.py'
METRICS=('disp_z_pct_hs','disp_3d_rms_m','roll_rms_deg','pitch_rms_deg','yaw_rms_deg','accel_bias_3d_rms_mps2','gyro_bias_3d_rms_radps')
WEIGHTS=dict(disp_z_pct_hs=1.5,disp_3d_rms_m=1.5,roll_rms_deg=.75,pitch_rms_deg=.75,yaw_rms_deg=.5,accel_bias_3d_rms_mps2=1.25,gyro_bias_3d_rms_radps=1.0)
BASE={'name':'baseline','env':{'TFG_TUNING':'adaptive','TFG_AW_COV_SYNC':'1'}}
def cfg(name,**env):
 e=dict(BASE['env']);e.update({k:str(v) for k,v in env.items()});return {'name':name,'env':e}
def screen():
 c=[BASE]
 axes={'TFG_TAU_COEFF':[.75,.9,1.1,1.25],'TFG_SIGMA_COEFF':[.6,.7,.9,1.0],
 'TFG_R_S_MSE_COEFF':[.035,.045,.065,.08],'TFG_R_S_X_FACTOR':[.8,1.0,1.3,1.55],
 'TFG_R_S_Y_FACTOR':[.8,1.0,1.3,1.55],'SF_SIGMA_G_SCALE':[.05,.15,.3,.6],
 'SF_SIGMA_M_SCALE':[2,3,5,7],'SF_GYRO_BIAS_RW_VAR':[3e-11,3e-10,1e-9],
 'TFG_ACC_BIAS_UNLOCK_SEC':[5,15,30],'TFG_ACC_BIAS_OU_TAU':[1500,3000,8000],
 'TFG_ADAPT_RS_MULT':[.6,1.0,2.5,4.0],'TFG_ACC_NOISE_FLOOR':[.08,.10,.15,.18]}
 for k,vals in axes.items():
  for v in vals:
   kw={k:v}
   if k=='TFG_ACC_BIAS_OU_TAU':kw['TFG_ACC_BIAS_OU_STD']=.025
   c.append(cfg(k.lower().replace('tfg_','').replace('sf_','')+'_'+str(v).replace('.','p'),**kw))
 for x,y in ((.9,.9),(1.0,1.0),(1.3,1.3),(1.5,1.2)):c.append(cfg(f'rsxy_{x}_{y}',TFG_R_S_X_FACTOR=x,TFG_R_S_Y_FACTOR=y))
 for g,m,b in ((.05,3,3e-11),(.1,4,1e-10),(.2,5,3e-10),(.3,6,1e-9)):c.append(cfg(f'attbias_{g}_{m}_{b}',SF_SIGMA_G_SCALE=g,SF_SIGMA_M_SCALE=m,SF_GYRO_BIAS_RW_VAR=b))
 return c
def run(stage,configs,seeds,out,jobs):
 d=out/stage;d.mkdir(parents=True,exist_ok=True);cp=d/'configs.json';cp.write_text(json.dumps(configs,indent=2)+'\n')
 subprocess.run([sys.executable,str(RUNNER),'--family','TFG','--configs',str(cp),'--records','0,1,2,3,4,5,6,7','--seeds',seeds,'--jobs',str(jobs),'--output-dir',str(d)],cwd=ROOT,check=True)
 return json.loads((d/'summary.json').read_text())
def objective(rows):
 by={r['config']:r for r in rows};base=by['baseline'];scores=[]
 for r in rows:
  parts={m:r[m]/max(base[m],1e-12) for m in METRICS}
  score=sum(WEIGHTS[m]*math.log(max(parts[m],1e-12)) for m in METRICS)/sum(WEIGHTS.values())
  penalty=.02*max(0,r['violations']-base['violations'])+.02*max(0,r['worst_ratio']-max(1,base['worst_ratio']))
  scores.append(dict(config=r['config'],score=score+penalty,raw_score=score,gate_penalty=penalty,ratios=parts,summary=r))
 return sorted(scores,key=lambda x:x['score'])
def refine(best):
 env={k:v for k,v in best['env'].items() if k not in BASE['env']};out=[BASE,{'name':'screen_winner','env':best['env']}]
 for k,v in list(env.items()):
  try:x=float(v)
  except ValueError:continue
  if x<=0:continue
  for f in (.8,.9,1.1,1.25):
   e=dict(best['env']);e[k]=str(x*f)
   if k=='TFG_ACC_BIAS_OU_TAU' and 'TFG_ACC_BIAS_OU_STD' not in e:e['TFG_ACC_BIAS_OU_STD']='0.025'
   out.append({'name':f'refine_{k.lower()}_{f}','env':e})
 for rs in (.85,1.15):
  for cov in (.75,1.25):
   e=dict(best['env'])
   if 'TFG_R_S_MSE_COEFF' in e:e['TFG_R_S_MSE_COEFF']=str(float(e['TFG_R_S_MSE_COEFF'])*rs)
   if 'SF_SIGMA_G_SCALE' in e:e['SF_SIGMA_G_SCALE']=str(float(e['SF_SIGMA_G_SCALE'])*cov)
   out.append({'name':f'cross_rs{rs}_cov{cov}','env':e})
 seen=set();uniq=[]
 for c in out:
  key=tuple(sorted(c['env'].items()))
  if key not in seen:seen.add(key);uniq.append(c)
 return uniq
def main():
 p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_rigorous_tuning');p.add_argument('--jobs',type=int,default=4);a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 s=screen();sr=run('screen',s,'default,11,23',out,a.jobs);rank=objective(sr);(out/'screen-ranking.json').write_text(json.dumps(rank,indent=2)+'\n')
 best=next(c for c in s if c['name']==rank[0]['config'])
 rc=refine(best);rr=run('refine',rc,'31,47,59',out,a.jobs);rrank=objective(rr);(out/'refine-ranking.json').write_text(json.dumps(rrank,indent=2)+'\n')
 final=next(c for c in rc if c['name']==rrank[0]['config']);final={'name':'candidate','env':final['env']}
 hold=run('holdout',[BASE,final],'101,1009,2027,3037',out,a.jobs)
 oucfg=out/'ou3-holdout-configs.json';oucfg.write_text(json.dumps([{'name':'shipping_ou3','env':{}}],indent=2)+'\n')
 subprocess.run([sys.executable,str(RUNNER),'--family','OU-III','--configs',str(oucfg),'--records','0,1,2,3,4,5,6,7','--seeds','101,1009,2027,3037','--jobs',str(a.jobs),'--output-dir',str(out/'ou3-holdout')],cwd=ROOT,check=True)
 ou=json.loads((out/'ou3-holdout/summary.json').read_text())[0]
 report={'objective_metrics':METRICS,'weights':WEIGHTS,'selection_rule':'minimum weighted log RMS ratio to stage baseline plus explicit gate-degradation penalty','screen_winner':rank[0],'refine_winner':rrank[0],'final_candidate':final,'holdout':hold,'ou3_holdout':ou}
 (out/'campaign-summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
