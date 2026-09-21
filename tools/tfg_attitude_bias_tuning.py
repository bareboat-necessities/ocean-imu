#!/usr/bin/env python3
"""Target TFG tilt/accelerometer-bias observability after the broad campaign."""
from __future__ import annotations
import argparse,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'tools/rao_parameter_tuning.py'
BASEENV={'TFG_TUNING':'adaptive','TFG_AW_COV_SYNC':'1','SF_GYRO_BIAS_RW_VAR':'3.75e-10'}
METRICS=('pitch_rms_deg','roll_rms_deg','accel_bias_3d_rms_mps2','disp_3d_rms_m','disp_z_pct_hs','yaw_rms_deg','gyro_bias_3d_rms_radps')
W=dict(pitch_rms_deg=2.0,roll_rms_deg=.75,accel_bias_3d_rms_mps2=2.0,disp_3d_rms_m=1.0,disp_z_pct_hs=1.0,yaw_rms_deg=.5,gyro_bias_3d_rms_radps=.5)
def arm(name,**kw):
 e=dict(BASEENV);e.update({k:str(v) for k,v in kw.items()});return {'name':name,'env':e}
def candidates():
 out=[arm('baseline')]
 axes={
 'TFG_ACC_BIAS_UNLOCK_SEC':[2,5,15,30,60],
 'TFG_HANDOFF_ACC_BIAS_STD':[.01,.02,.05,.08,.12],
 'TFG_HANDOFF_TILT_SIGMA_RAD':[.012,.02,.05,.08,.12],
 'TFG_RACC_WARMUP_STD':[.15,.25,.35,.7,1.0],
 'TFG_PROXY_TWO_KP':[.1,.15,.3,.4,.6],
 'TFG_PROXY_TWO_KI':[.005,.01,.04,.08],
 'TFG_PROXY_GRAVITY_LPF_SEC':[6,9,18,24],
 'TFG_ACC_BIAS_OU_TAU':[1000,2000,3500,7000,12000],
 }
 for k,vs in axes.items():
  for v in vs:
   kw={k:v}
   if k=='TFG_ACC_BIAS_OU_TAU':kw['TFG_ACC_BIAS_OU_STD']=.025
   out.append(arm(k.lower().replace('tfg_','')+'_'+str(v).replace('.','p'),**kw))
 # Joint covariance/observability arms: tilt confidence x bias freedom.
 for tilt in (.015,.025,.05,.08):
  for bs in (.015,.03,.06,.1):
   out.append(arm(f'joint_tilt{tilt}_bias{bs}',TFG_HANDOFF_TILT_SIGMA_RAD=tilt,TFG_HANDOFF_ACC_BIAS_STD=bs))
 for unlock,tau in ((3,2000),(5,3500),(15,3500),(30,7000),(60,12000)):
  out.append(arm(f'biasdyn_u{unlock}_t{tau}',TFG_ACC_BIAS_UNLOCK_SEC=unlock,TFG_ACC_BIAS_OU_TAU=tau,TFG_ACC_BIAS_OU_STD=.025))
 return out
def execute(tag,cfgs,seeds,out,jobs):
 d=out/tag;d.mkdir(parents=True,exist_ok=True);p=d/'configs.json';p.write_text(json.dumps(cfgs,indent=2)+'\n')
 subprocess.run([sys.executable,str(RUN),'--family','TFG','--configs',str(p),'--records','0,1,2,3,4,5,6,7','--seeds',seeds,'--jobs',str(jobs),'--output-dir',str(d)],cwd=ROOT,check=True)
 return json.loads((d/'summary.json').read_text())
def rank(rows):
 b=next(r for r in rows if r['config']=='baseline');z=[]
 for r in rows:
  ratios={m:r[m]/max(b[m],1e-12) for m in METRICS}
  score=sum(W[m]*math.log(max(ratios[m],1e-12)) for m in METRICS)/sum(W.values())
  penalty=.03*max(0,r['violations']-b['violations'])+.03*max(0,r['worst_ratio']-max(1,b['worst_ratio']))
  z.append({'config':r['config'],'score':score+penalty,'ratios':ratios,'summary':r})
 return sorted(z,key=lambda x:x['score'])
def refine(best):
 out=[arm('baseline'),{'name':'screen_winner','env':best['env']}]
 # Local multiplicative refinement of every selected non-baseline knob.
 for k,v in best['env'].items():
  if k in BASEENV:continue
  try:x=float(v)
  except:continue
  for f in (.7,.85,1.15,1.4):
   e=dict(best['env']);e[k]=str(x*f)
   if k=='TFG_ACC_BIAS_OU_TAU':e.setdefault('TFG_ACC_BIAS_OU_STD','.025')
   out.append({'name':f'local_{k.lower()}_{f}','env':e})
 # Deliberately combine best with likely covariance separators.
 for tilt in (.015,.03,.06):
  for bs in (.015,.04,.08):
   e=dict(best['env'],TFG_HANDOFF_TILT_SIGMA_RAD=str(tilt),TFG_HANDOFF_ACC_BIAS_STD=str(bs))
   out.append({'name':f'cross_t{tilt}_b{bs}','env':e})
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_attitude_bias_tuning');p.add_argument('--jobs',type=int,default=4);a=p.parse_args();o=a.output_dir.resolve();o.mkdir(parents=True,exist_ok=True)
 cs=candidates();r=execute('screen',cs,'default,11,23',o,a.jobs);q=rank(r);(o/'screen-ranking.json').write_text(json.dumps(q,indent=2)+'\n');best=next(c for c in cs if c['name']==q[0]['config'])
 cr=refine(best);rr=execute('refine',cr,'31,47,59',o,a.jobs);qr=rank(rr);(o/'refine-ranking.json').write_text(json.dumps(qr,indent=2)+'\n');winner=next(c for c in cr if c['name']==qr[0]['config']);winner={'name':'candidate','env':winner['env']}
 h=execute('holdout',[arm('baseline'),winner],'101,1009,2027,3037',o,a.jobs)
 summary={'focus':'physical tilt/accelerometer-bias observability and covariance; no gate changes','metrics':METRICS,'weights':W,'screen_winner':q[0],'refine_winner':qr[0],'candidate':winner,'holdout':h}
 (o/'campaign-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
