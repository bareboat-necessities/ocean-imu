#!/usr/bin/env python3
import argparse, csv, math, os, random, re, subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports' / 'results'

GRAVITY_RANGES = {
    'OU_II': {
        'SF_BOOT_TILT_ACC_TAU': (1.5, 2.8), 'SF_BOOT_GRAV_SLOW_TAU': (3.5, 6.0), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.12, 0.20),
        'SF_BOOT_GRAV_HOLD_SEC': (0.6, 1.5), 'SF_BOOT_GRAV_MIN_SEC': (4.0, 8.0), 'SF_BOOT_GRAV_NORM_FRAC': (0.30, 0.49),
        'SF_ONLINE_TUNE_WARMUP_SEC': (5.0, 10.0),
    },
    'OU_III': {
        'SF_BOOT_TILT_ACC_TAU': (1.5, 3.2), 'SF_BOOT_GRAV_SLOW_TAU': (3.5, 7.0), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.10, 0.20),
        'SF_BOOT_GRAV_HOLD_SEC': (0.6, 2.0), 'SF_BOOT_GRAV_MIN_SEC': (4.0, 9.0), 'SF_BOOT_GRAV_NORM_FRAC': (0.25, 0.49),
        'SF_ONLINE_TUNE_WARMUP_SEC': (5.0, 12.0),
    }
}
OU_COMMON = {
    'OU_ACC_NOISE_FLOOR_SIGMA': (0.10, 0.18), 'OU_SIGMA_COEFF': (0.70, 0.95), 'OU_TAU_COEFF': (1.3, 1.8),
    'OU_ADAPT_TAU_SEC': (1.5, 3.5), 'OU_ADAPT_EVERY_SECS': (0.08, 0.20)
}
OU_II_ONLY = {'OU_P_FACTOR': (1.2, 1.8), 'OU_R_P0_XY_FACTOR': (0.20, 0.40), 'OU_R_P0_COEFF': (1.2, 2.2), 'OU_R_V0_COEFF': (1.1, 2.0)}

RX={'ds':re.compile(r"Processing\s+(.+?\.csv)"),'ang':re.compile(r"Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw=([0-9.eE+-]+)"),'xyz':re.compile(r"XYZ RMS \(m\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+)"),'r3d':re.compile(r"3D RMS \(m\):\s*([0-9.eE+-]+)"),'accb':re.compile(r"Bias error RMS \(acc, m/s\^2\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+) \|3D\|=([0-9.eE+-]+)"),'gate':re.compile(r"QUALITY_GATE:\s+PASS=(\d)\s+REASON=(.*)")}

def lhs(r,n,rng):
 a=[rng[0]+(i+r.random())/n*(rng[1]-rng[0]) for i in range(n)]; r.shuffle(a); return a

def minimal_candidates():
 return [(f'minimal_{i:03d}',p) for i,p in enumerate([
    {'SF_BOOT_GRAV_HOLD_SEC':1.0},{'SF_BOOT_GRAV_MIN_SEC':6.0},{'SF_BOOT_GRAV_NORM_FRAC':0.35},{'SF_BOOT_GRAV_SLOW_TAU':5.0},{'SF_BOOT_TILT_ACC_TAU':2.0},
    {'SF_ONLINE_TUNE_WARMUP_SEC':8.0},{'SF_BOOT_GRAV_HOLD_SEC':1.0,'SF_BOOT_GRAV_MIN_SEC':6.0},{'SF_BOOT_GRAV_NORM_FRAC':0.35,'SF_BOOT_GRAV_SLOW_TAU':5.0},
    {'SF_BOOT_TILT_ACC_TAU':2.0,'SF_BOOT_GRAV_HOLD_SEC':1.0},{'SF_ONLINE_TUNE_WARMUP_SEC':8.0,'SF_BOOT_GRAV_MIN_SEC':6.0}], start=1)]

def sample_params(space, n, r):
 S={k:lhs(r,n,v) for k,v in space.items()}; out=[]
 for i in range(n):
  p={k:S[k][i] for k in S}
  margin=2.0+r.random()*6.0; p['SF_BOOT_GRAV_TIMEOUT_SEC']=p['SF_BOOT_GRAV_MIN_SEC']+p['SF_BOOT_GRAV_HOLD_SEC']+margin
  out.append(p)
 return out

def parse(text):
 ds=RX['ds'].findall(text); ang=[m.groups() for m in RX['ang'].finditer(text)]; xyz=[m.groups() for m in RX['xyz'].finditer(text)]; r3d=[m.group(1) for m in RX['r3d'].finditer(text)]; accb=[m.groups() for m in RX['accb'].finditer(text)]; gates=[m.groups() for m in RX['gate'].finditer(text)]
 n=max(len(ds),len(ang),len(xyz),len(accb),1); out=[]
 for i in range(n):
  g=lambda arr,j,idx=float('nan'): float(arr[j][idx]) if j<len(arr) else float('nan')
  gp=(gates[i][0]=='1') if i<len(gates) else False; reason=(gates[i][1].strip() if i<len(gates) else 'missing_gate_status')
  out.append(dict(wave_dataset=ds[i] if i<len(ds) else f'dataset_{i}',roll_rms=g(ang,i,0),pitch_rms=g(ang,i,1),yaw_rms=g(ang,i,2),x_rms=g(xyz,i,0),y_rms=g(xyz,i,1),z_rms=g(xyz,i,2),rms_3d=float(r3d[i]) if i<len(r3d) else float('nan'),acc_bias_rms_3d=g(accb,i,3),quality_gate_pass=gp,fail_reason=reason))
 return out

def run_candidate(fam,cid,p,tier,seed,collect):
 tdir=ROOT/'tests'/('kalman_ou_ii' if fam=='OU_II' else 'kalman_ou_iii'); bin_name='./kalman_ou_ii-sim' if fam=='OU_II' else './kalman_ou_iii-sim'
 env=os.environ.copy(); env['W3D_SEED']=str(seed); env['W3D_TIER']=tier
 if collect: env['W3D_COLLECT_ALL_GATES']='1'
 env.update({k:(f'{v:.6g}' if isinstance(v,float) else str(v)) for k,v in p.items()})
 pr=subprocess.run([bin_name,'--nomag'],cwd=tdir,text=True,capture_output=True,env=env)
 rows=parse(pr.stdout+'\n'+pr.stderr)
 for r in rows:
  r.update({'family':fam,'candidate':cid,'seed':seed,'tier':tier,'returncode':pr.returncode,'roll_pitch_rms_norm':math.hypot(r['roll_rms'],r['pitch_rms']),'xy_rms':math.hypot(r['x_rms'],r['y_rms'])})
  r.update(p)
 return rows

def eval_set(fam, cand_list, tier='quick', collect=True, seed=42):
 rows=[]; total=len(cand_list); report_every=max(1,total//10)
 for i,(cid,p) in enumerate(cand_list,1):
  rows.extend(run_candidate(fam,cid,p,tier,seed,collect))
  if i==1 or i==total or i%report_every==0:
   print(f'[{fam}] progress: {i}/{total} candidates complete')
 return rows

def percentile(vals,p):
 if not vals:return float('inf')
 s=sorted(vals); i=(len(s)-1)*p; lo=int(i); hi=min(len(s)-1,lo+1); f=i-lo; return s[lo]*(1-f)+s[hi]*f

def score_rows(rows):
 by=defaultdict(list)
 for r in rows: by[(r['family'],r['candidate'],r['tier'])].append(r)
 base={(r['family'],r['wave_dataset'],r['seed'],r['tier']):r for r in rows if r['candidate']=='baseline'}
 for _,it in by.items():
  rp=[]; z=[]; r3=[]; yaw=[]; acc=[]; any_fail=False
  for r in it:
   b=base.get((r['family'],r['wave_dataset'],r['seed'],r['tier']),r)
   rf=lambda x: x if math.isfinite(x) and abs(x)>1e-9 else 1e-9
   r['baseline_roll_pitch_rms_norm']=b['roll_pitch_rms_norm']; r['baseline_z_rms']=b['z_rms']; r['baseline_rms_3d']=b['rms_3d']; r['baseline_yaw_rms']=b['yaw_rms']
   r['roll_pitch_ratio']=rf(r['roll_pitch_rms_norm'])/rf(b['roll_pitch_rms_norm']); r['z_ratio']=rf(r['z_rms'])/rf(b['z_rms']); r['rms3d_ratio']=rf(r['rms_3d'])/rf(b['rms_3d']); r['yaw_ratio']=rf(r['yaw_rms'])/rf(b['yaw_rms']); r['acc_bias_ratio']=rf(r['acc_bias_rms_3d'])/rf(b['acc_bias_rms_3d'])
   rp.append(r['roll_pitch_ratio']); z.append(r['z_ratio']); r3.append(r['rms3d_ratio']); yaw.append(r['yaw_ratio']); acc.append(r['acc_bias_ratio']); any_fail = any_fail or (not r['quality_gate_pass'])
  s=4.0*percentile(rp,0.75)+2.0*max(rp)+2.0*max(0,max(z)-1.02)+2.0*max(0,max(r3)-1.02)+0.5*percentile(yaw,0.75)+0.5*percentile(acc,0.75)+(100.0 if any_fail else 0.0)
  for r in it: r['score']=s; r['survived_tier']=int((not any_fail) and max(z)<=1.05 and max(r3)<=1.05)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--mode',choices=['gravity','ou','full'],default='gravity'); ap.add_argument('--family',choices=['OU_II','OU_III','both'],default='both'); ap.add_argument('--samples',type=int,default=100); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--tier',choices=['quick','all','final'],default='quick'); ap.add_argument('--top-k',type=int,default=8); ap.add_argument('--collect-all-gates',action='store_true',default=False)
 a=ap.parse_args(); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_ii'),'build'],check=True); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_iii'),'build'],check=True)
 fams=['OU_II','OU_III'] if a.family=='both' else [a.family]; rows=[]
 for fam in fams:
  rng=random.Random(a.seed+(0 if fam=='OU_II' else 1000000)); grav_mc=[(f'grav_mc_{i:04d}',p) for i,p in enumerate(sample_params(GRAVITY_RANGES[fam],a.samples,rng),1)]
  candidates=[('baseline',{})]+minimal_candidates()+grav_mc
  if a.mode in ('ou','full'):
   ou_space=dict(OU_COMMON); 
   if fam=='OU_II': ou_space.update(OU_II_ONLY)
   candidates += [(f'ou_mc_{i:04d}',p) for i,p in enumerate(sample_params(ou_space,a.samples,rng),1)]
  rows.extend(eval_set(fam,candidates,tier=a.tier,collect=(a.collect_all_gates or a.tier!='final'),seed=a.seed))
 score_rows(rows); REPORTS.mkdir(parents=True,exist_ok=True)
 out=REPORTS/f'ou_sweep_{a.mode}_{a.tier}_seed{a.seed}_n{a.samples}.csv'
 fields=sorted({k for r in rows for k in r.keys()})
 with out.open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
 print(out)

if __name__=='__main__': main()
