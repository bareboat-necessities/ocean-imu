#!/usr/bin/env python3
import argparse, csv, math, os, random, re, subprocess
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; REPORTS=ROOT/'reports'/'results'

COMMON={"SF_BOOT_TILT_ACC_TAU":(1.5,5.0),"SF_BOOT_GRAV_SLOW_TAU":(3.5,10.0),"SF_BOOT_GRAV_ALIGN_MAX_SIN":(0.04,0.12),"SF_BOOT_GRAV_HOLD_SEC":(1.0,5.0),"SF_BOOT_GRAV_MIN_SEC":(5.0,16.0),"SF_BOOT_GRAV_TIMEOUT_SEC":(10.0,30.0),"SF_BOOT_GRAV_NORM_FRAC":(0.15,0.40),"SF_ONLINE_TUNE_WARMUP_SEC":(8.0,20.0),"OU_ACC_NOISE_FLOOR_SIGMA":(0.10,0.25),"OU_SIGMA_COEFF":(0.60,1.00),"OU_TAU_COEFF":(1.2,2.2),"OU_ADAPT_TAU_SEC":(1.5,5.0),"OU_ADAPT_EVERY_SECS":(0.05,0.25)}
OU2={"OU_P_FACTOR":(1.0,2.2),"OU_R_P0_XY_FACTOR":(0.10,0.60),"OU_R_P0_COEFF":(1.0,3.0),"OU_R_V0_COEFF":(1.0,3.0)}
MAG={"SF_MAG_DELAY_SEC":(5.0,25.0),"SF_MAG_GRAV_ALIGN_MAX_SIN":(0.04,0.16),"SF_MAG_GRAV_ALIGN_HOLD_SEC":(0.2,4.0),"SF_MAG_GRAV_ALIGN_LPF_TAU":(0.2,2.0),"SF_MAG_TILT_FALLBACK_SEC":(3.0,40.0),"SF_MAG_EXTREME_GYRO_DPS":(35.0,130.0),"SF_MAG_MIN_SAMPLES":(24,800),"SF_MAG_INIT_MIN_MAG_NORM":(1e-4,5e-3)}

RX={
'ds':re.compile(r"Processing\s+(.+?\.csv)"),
'ang':re.compile(r"Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw=([0-9.eE+-]+)"),
'xyz':re.compile(r"XYZ RMS \(m\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+)"),
'r3d':re.compile(r"3D RMS \(m\):\s*([0-9.eE+-]+)"),
'accb':re.compile(r"Bias error RMS \(acc, m/s\^2\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+) \|3D\|=([0-9.eE+-]+)")}

def lhs(r,n,rng):
 a=[rng[0]+(i+r.random())/n*(rng[1]-rng[0]) for i in range(n)]; r.shuffle(a); return a

def cands(fam,samples,seed,with_mag):
 r=random.Random(seed+(0 if fam=='OU_II' else 1_000_000)); out=[('baseline',{})]
 explicit={"SF_BOOT_TILT_ACC_TAU":2.5,"SF_BOOT_GRAV_SLOW_TAU":6.0,"SF_BOOT_GRAV_ALIGN_MAX_SIN":0.070,"SF_BOOT_GRAV_HOLD_SEC":2.0,"SF_BOOT_GRAV_MIN_SEC":8.0,"SF_BOOT_GRAV_TIMEOUT_SEC":15.0,"SF_BOOT_GRAV_NORM_FRAC":0.25,"SF_ONLINE_TUNE_WARMUP_SEC":10.0}
 for rr in (0.8,1.2,1.6): out.append((f'explicit_racc_{str(rr).replace('.','_')}',{**explicit,'SF_RACC_WARMUP_STD':rr}))
 for rr in (0.8,1.2,1.6): out.append((f'racc_{str(rr).replace('.','_')}',{'SF_RACC_WARMUP_STD':rr}))
 space=dict(COMMON); 
 if fam=='OU_II': space.update(OU2)
 if with_mag: space.update(MAG)
 S={k:lhs(r,samples,v) for k,v in space.items()}
 for i in range(samples):
  p={k:S[k][i] for k in S}; p['SF_RACC_WARMUP_STD']=(0.8,1.2,1.6)[i%3];
  if 'SF_MAG_MIN_SAMPLES' in p: p['SF_MAG_MIN_SAMPLES']=int(round(p['SF_MAG_MIN_SAMPLES']))
  out.append((f'mc_{i:04d}',p))
 return out

def parse(text):
 ds=RX['ds'].findall(text); ang=[m.groups() for m in RX['ang'].finditer(text)]; xyz=[m.groups() for m in RX['xyz'].finditer(text)]; r3d=[m.group(1) for m in RX['r3d'].finditer(text)]; accb=[m.groups() for m in RX['accb'].finditer(text)]
 n=max(len(ds),len(ang),len(xyz),len(accb),1); out=[]
 for i in range(n):
  g=lambda arr,j,idx=float('nan'): float(arr[j][idx]) if j<len(arr) else float('nan')
  out.append(dict(wave_dataset=ds[i] if i<len(ds) else f'dataset_{i}',roll_rms=g(ang,i,0),pitch_rms=g(ang,i,1),yaw_rms=g(ang,i,2),x_rms=g(xyz,i,0),y_rms=g(xyz,i,1),z_rms=g(xyz,i,2),rms_3d=float(r3d[i]) if i<len(r3d) else float('nan'),acc_bias_rms_x=g(accb,i,0),acc_bias_rms_y=g(accb,i,1),acc_bias_rms_z=g(accb,i,2),acc_bias_rms_3d=g(accb,i,3)))
 return out

def run_family(fam,args):
 tdir=ROOT/'tests'/('kalman_ou_ii' if fam=='OU_II' else 'kalman_ou_iii'); bin_name='./kalman_ou_ii-sim' if fam=='OU_II' else './kalman_ou_iii-sim'; rows=[]
 cc=cands(fam,args.samples,args.seed,args.with_mag_sweep); total=len(cc); interval=max(1,total//20)
 for idx,(cid,p) in enumerate(cc,start=1):
  env=os.environ.copy(); env.update({k:f'{v:.6g}' if isinstance(v,float) else str(v) for k,v in p.items()})
  pr=subprocess.run([bin_name],cwd=tdir,text=True,capture_output=True,env=env)
  for m in parse(pr.stdout+'\n'+pr.stderr):
   m.update(dict(family=fam,candidate=cid,seed=args.seed,returncode=pr.returncode,quality_gate_pass=(pr.returncode==0),fail_reason='' if pr.returncode==0 else 'nonzero_return_or_gate_fail',roll_pitch_rms_norm=math.hypot(m['roll_rms'],m['pitch_rms']) if math.isfinite(m['roll_rms']) and math.isfinite(m['pitch_rms']) else float('nan'),xy_rms=math.hypot(m['x_rms'],m['y_rms']) if math.isfinite(m['x_rms']) and math.isfinite(m['y_rms']) else float('nan')))
   for k,v in p.items(): m[k]=v
   rows.append(m)
  if idx==1 or idx==total or (idx%interval)==0:
   pct=100.0*idx/total if total else 100.0
   print(f'[{fam}] progress {idx}/{total} ({pct:.1f}%)',flush=True)
 return rows

def score(rows):
 base={(r['family'],r['wave_dataset'],r['seed']):r for r in rows if r['candidate']=='baseline'}
 for r in rows:
  b=base.get((r['family'],r['wave_dataset'],r['seed']),r)
  def fv(x,d=100.0):
   try:v=float(x); return v if math.isfinite(v) else d
   except:return d
  acc=fv(r['acc_bias_rms_3d']); rp=fv(r['roll_pitch_rms_norm']); yaw=fv(r['yaw_rms']); z=fv(r['z_rms']); xy=fv(r['xy_rms']); bz=fv(b['z_rms']); bxy=fv(b['xy_rms'])
  r['score']=8*acc+3*rp+1*yaw+2*max(0,z-bz)+2*max(0,xy-bxy)+(0 if r['quality_gate_pass'] else 1e6)+(0 if int(r['returncode'])==0 else 1e6)

def report(rows,args):
 REPORTS.mkdir(parents=True,exist_ok=True)
 out=REPORTS/f'ou_sweep_seed{args.seed}_n{args.samples}_mag{int(args.with_mag_sweep)}.csv'
 first_fields=list(rows[0].keys())
 extra_fields=sorted({k for r in rows for k in r.keys()}-set(first_fields))
 all_fields=first_fields+extra_fields
 with out.open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=all_fields); w.writeheader(); w.writerows(rows)
 for fam in sorted({r['family'] for r in rows}):
  fr=[r for r in rows if r['family']==fam]; by=defaultdict(list)
  for r in fr: by[r['candidate']].append(r)
  agg=[]
  for c,it in by.items():
   if not all(bool(x.get('quality_gate_pass')) and int(x.get('returncode',1))==0 for x in it): continue
   agg.append((sum(float(x['score']) for x in it)/len(it),c,it))
  if not agg:
   raise RuntimeError(f'No quality-gate-passing candidates found for {fam}')
  best=min(agg,key=lambda x:x[0])
  txt=REPORTS/f'{fam.lower()}_best_report.txt'
  with txt.open('w') as f:
   f.write(f'family: {fam}\nbest_candidate: {best[1]}\nmean_score: {best[0]:.6f}\nquality_rows: {sum(1 for x in best[2] if x["quality_gate_pass"])} / {len(best[2])}\n')
   f.write('priority: accel_bias_rms_3d > attitude(r/p norm) > yaw\n')
   f.write('\n')
   param_keys=sorted({k for k in best[2][0].keys() if k.startswith('SF_') or k.startswith('OU_')})
   f.write('best_found_parameters:\n')
   for k in param_keys: f.write(f'  {k}: {best[2][0][k]}\n')
   f.write('\n')
   f.write('per_wave_rms:\n')
   for r in sorted(best[2],key=lambda x:x['wave_dataset']):
    f.write(f"  wave_dataset: {r['wave_dataset']}\n")
    f.write(f"    xyz_rms_m: X={r['x_rms']}, Y={r['y_rms']}, Z={r['z_rms']}\n")
    f.write(f"    rms_3d_m: {r['rms_3d']}\n")
    f.write(f"    angles_rms_deg: Roll={r['roll_rms']}, Pitch={r['pitch_rms']}, Yaw={r['yaw_rms']}\n")
    f.write(f"    acc_bias_rms_mps2: X={r['acc_bias_rms_x']}, Y={r['acc_bias_rms_y']}, Z={r['acc_bias_rms_z']}, |3D|={r['acc_bias_rms_3d']}\n")
    f.write(f"    quality_gate_pass: {r['quality_gate_pass']}\n")
 return out

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=100); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--family',choices=['OU_II','OU_III','both'],default='both'); ap.add_argument('--with-mag-sweep',action='store_true',default=False)
 a=ap.parse_args(); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_ii'),'build'],check=True); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_iii'),'build'],check=True)
 fams=['OU_II','OU_III'] if a.family=='both' else [a.family]; rows=[]
 for fam in fams: rows.extend(run_family(fam,a))
 score(rows); out=report(rows,a); print(out)

if __name__=='__main__': main()
