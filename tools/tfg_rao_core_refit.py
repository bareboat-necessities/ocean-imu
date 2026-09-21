#!/usr/bin/env python3
"""RAO-only refit of TFG c_tau, c_sigma and horizontal S regularization.\n

This supersedes historical comments fitted before the vessel-RAO dataset.
It is deliberately narrow: shipping gyro-bias tuning is held fixed, no other
filter coefficients move, and the old long generic campaigns are not rerun.
"""
from __future__ import annotations
import argparse,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'tools/rao_parameter_tuning.py'
BASEENV={'TFG_TUNING':'adaptive','TFG_AW_COV_SYNC':'1'}
METRICS=('disp_z_pct_hs','disp_3d_rms_m','roll_rms_deg','pitch_rms_deg','yaw_rms_deg','accel_bias_3d_rms_mps2','gyro_bias_3d_rms_radps')
W=dict(disp_z_pct_hs=2.,disp_3d_rms_m=2.,roll_rms_deg=.5,pitch_rms_deg=.5,yaw_rms_deg=.25,accel_bias_3d_rms_mps2=.75,gyro_bias_3d_rms_radps=.5)
def arm(name,tau=1.,sigma=.8,x=1.15,y=1.15):
 return {'name':name,'env':dict(BASEENV,TFG_TAU_COEFF=str(tau),TFG_SIGMA_COEFF=str(sigma),TFG_R_S_X_FACTOR=str(x),TFG_R_S_Y_FACTOR=str(y))}
def run(tag,cfgs,seeds,out,jobs):
 d=out/tag;d.mkdir(parents=True,exist_ok=True);p=d/'configs.json';p.write_text(json.dumps(cfgs,indent=2)+'\n')
 subprocess.run([sys.executable,str(RUN),'--family','TFG','--configs',str(p),'--records','0,1,2,3,4,5,6,7','--seeds',seeds,'--jobs',str(jobs),'--output-dir',str(d)],cwd=ROOT,check=True)
 return json.loads((d/'summary.json').read_text())
def rank(rows):
 b=next(r for r in rows if r['config']=='baseline');ans=[]
 for r in rows:
  ratios={m:r[m]/max(b[m],1e-12) for m in METRICS}
  score=sum(W[m]*math.log(max(ratios[m],1e-12)) for m in METRICS)/sum(W.values())
  penalty=.025*max(0,r['violations']-b['violations'])+.025*max(0,r['worst_ratio']-max(1,b['worst_ratio']))
  ans.append(dict(config=r['config'],score=score+penalty,raw_score=score,penalty=penalty,ratios=ratios,summary=r))
 return sorted(ans,key=lambda z:z['score'])
def screen():
 out=[arm('baseline')]
 # Wide one-dimensional RAO refit: old defaults are not assumed near optimum.
 for v in (.45,.6,.75,.9,1.1,1.3,1.55,1.85):out.append(arm(f'tau_{v}',tau=v))
 for v in (.35,.5,.65,.95,1.15,1.4):out.append(arm(f'sigma_{v}',sigma=v))
 for v in (.45,.65,.85,1.0,1.3,1.55,1.85,2.2):out.append(arm(f'rsxy_{v}',x=v,y=v))
 # Coarse joint surface catches interactions hidden by 1-D sweeps.
 for t in (.65,.9,1.15,1.45):
  for s in (.5,.75,1.0,1.25):
   for xy in (.65,1.0,1.35,1.75):out.append(arm(f'joint_t{t}_s{s}_xy{xy}',t,s,xy,xy))
 return out
def refine(best):
 e=best['env'];t=float(e.get('TFG_TAU_COEFF',1));s=float(e.get('TFG_SIGMA_COEFF',.8));x=float(e.get('TFG_R_S_X_FACTOR',1.15));y=float(e.get('TFG_R_S_Y_FACTOR',1.15))
 vals=lambda z:[z*.78,z*.9,z,z*1.1,z*1.22]
 out=[arm('baseline')]
 for tt in vals(t):
  for ss in vals(s):
   for xx in vals(x):
    out.append(arm(f'r_t{tt:.4g}_s{ss:.4g}_x{xx:.4g}',tt,ss,xx,xx))
 # Explicit anisotropy around the symmetric winner.
 for fx,fy in ((.8,1),(1,.8),(1,1.2),(1.2,1),(.9,1.1),(1.1,.9)):
  out.append(arm(f'aniso_{fx}_{fy}',t,s,x*fx,y*fy))
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_rao_core_refit');p.add_argument('--jobs',type=int,default=4);a=p.parse_args();o=a.output_dir.resolve();o.mkdir(parents=True,exist_ok=True)
 cs=screen();q=rank(run('screen',cs,'default,11,23',o,a.jobs));(o/'screen-ranking.json').write_text(json.dumps(q,indent=2)+'\n');best=next(c for c in cs if c['name']==q[0]['config'])
 cr=refine(best);qr=rank(run('refine',cr,'31,47,59',o,a.jobs));(o/'refine-ranking.json').write_text(json.dumps(qr,indent=2)+'\n');win=next(c for c in cr if c['name']==qr[0]['config']);win={'name':'candidate','env':win['env']}
 hold=run('holdout',[arm('baseline'),win],'101,1009,2027,3037',o,a.jobs)
 (o/'campaign-summary.json').write_text(json.dumps({'purpose':'RAO refit; historical pre-RAO comments are not evidence','screen_winner':q[0],'refine_winner':qr[0],'candidate':win,'holdout':hold},indent=2)+'\n')
if __name__=='__main__':main()
