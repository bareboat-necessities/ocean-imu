#!/usr/bin/env python3
import argparse, csv, math, os, random, re, subprocess, time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_BOOT_TIMEOUT_SEC = 15.0
MIN_MARGIN_SEC = 1.0
MAX_MARGIN_SEC = 4.0
REJECTED_KEYS = {'OU_TAU_COEFF', 'OU_SIGMA_COEFF'}
REJECTED_PREFIXES = ('SF_MAG_',)

GRAVITY_RANGES = {
    'OU_II': {
        'SF_BOOT_TILT_ACC_TAU': (1.0, 4.0), 'SF_BOOT_GRAV_SLOW_TAU': (2.5, 8.0), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.08, 0.25),
        'SF_BOOT_GRAV_HOLD_SEC': (0.4, 2.2), 'SF_BOOT_GRAV_MIN_SEC': (3.0, 10.0), 'SF_BOOT_GRAV_NORM_FRAC': (0.18, 0.55),
        'SF_ONLINE_TUNE_WARMUP_SEC': (3.0, 15.0),
    },
    'OU_III': {
        'SF_BOOT_TILT_ACC_TAU': (1.0, 4.0), 'SF_BOOT_GRAV_SLOW_TAU': (2.5, 8.0), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.08, 0.25),
        'SF_BOOT_GRAV_HOLD_SEC': (0.4, 2.2), 'SF_BOOT_GRAV_MIN_SEC': (3.0, 10.0), 'SF_BOOT_GRAV_NORM_FRAC': (0.18, 0.55),
        'SF_ONLINE_TUNE_WARMUP_SEC': (3.0, 15.0),
    }
}
OU_COMMON = {
    'OU_ACC_NOISE_FLOOR_SIGMA': (0.06, 0.25),
    'OU_ADAPT_TAU_SEC': (0.8, 6.0), 'OU_ADAPT_EVERY_SECS': (0.04, 0.35)
}
OU_II_ONLY = {'OU_P_FACTOR': (0.8, 2.5), 'OU_R_P0_XY_FACTOR': (0.10, 0.70), 'OU_R_P0_COEFF': (0.6, 3.5), 'OU_R_V0_COEFF': (0.6, 3.0)}

RX={'ds':re.compile(r"Processing\s+(.+?\.csv)"),'ang':re.compile(r"Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw=([0-9.eE+-]+)"),'xyz':re.compile(r"XYZ RMS \(m\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+)"),'r3d':re.compile(r"3D RMS \(m\):\s*([0-9.eE+-]+)"),'accb':re.compile(r"Bias error RMS \(acc, m/s\^2\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+) \|3D\|=([0-9.eE+-]+)"),'gate':re.compile(r"QUALITY_GATE:\s+PASS=(\d)\s+REASON=(.*)")}

def lhs(r,n,rng):
 a=[rng[0]+(i+r.random())/n*(rng[1]-rng[0]) for i in range(n)]; r.shuffle(a); return a

def minimal_candidates():
 return [(f'minimal_{i:03d}',p) for i,p in enumerate([
    {'SF_BOOT_GRAV_HOLD_SEC':1.0},{'SF_BOOT_GRAV_MIN_SEC':6.0},{'SF_BOOT_GRAV_NORM_FRAC':0.35},{'SF_BOOT_GRAV_SLOW_TAU':5.0},{'SF_BOOT_TILT_ACC_TAU':2.0},
    {'SF_ONLINE_TUNE_WARMUP_SEC':8.0},{'SF_BOOT_GRAV_HOLD_SEC':1.0,'SF_BOOT_GRAV_MIN_SEC':6.0},{'SF_BOOT_GRAV_NORM_FRAC':0.35,'SF_BOOT_GRAV_SLOW_TAU':5.0},
    {'SF_BOOT_TILT_ACC_TAU':2.0,'SF_BOOT_GRAV_HOLD_SEC':1.0},{'SF_ONLINE_TUNE_WARMUP_SEC':8.0,'SF_BOOT_GRAV_MIN_SEC':6.0}], start=1)]

def explicit_practical_candidates():
 base = {
   'SF_BOOT_TILT_ACC_TAU':2.5, 'SF_BOOT_GRAV_SLOW_TAU':5.5, 'SF_BOOT_GRAV_ALIGN_MAX_SIN':0.12,
   'SF_BOOT_GRAV_HOLD_SEC':1.2, 'SF_BOOT_GRAV_MIN_SEC':6.0, 'SF_BOOT_GRAV_NORM_FRAC':0.35,
   'SF_ONLINE_TUNE_WARMUP_SEC':8.0
 }
 out = []
 for i, timeout in enumerate((8.0, 10.0, 12.0, 15.0), 1):
  p = dict(base)
  p['SF_BOOT_GRAV_TIMEOUT_SEC'] = timeout
  out.append((f'practical_{i:03d}', p))
 return out

def apply_timeout_constraint(p, r):
 min_sec = p['SF_BOOT_GRAV_MIN_SEC']
 hold_sec = p['SF_BOOT_GRAV_HOLD_SEC']
 margin_max = min(MAX_MARGIN_SEC, MAX_BOOT_TIMEOUT_SEC - min_sec - hold_sec)
 if margin_max < MIN_MARGIN_SEC:
  return False
 margin = MIN_MARGIN_SEC + r.random() * (margin_max - MIN_MARGIN_SEC)
 p['SF_BOOT_GRAV_TIMEOUT_SEC'] = min_sec + hold_sec + margin
 return p['SF_BOOT_GRAV_TIMEOUT_SEC'] <= MAX_BOOT_TIMEOUT_SEC

def sample_params(space, n, r, gravity=False):
 S={k:lhs(r,n,v) for k,v in space.items()}; out=[]
 i = 0
 while i < n:
  p={k:S[k][i] for k in S}
  if gravity and not apply_timeout_constraint(p, r):
   continue
  out.append(p)
  i += 1
 return out

def validate_candidate(p):
 bad = []
 timeout = p.get('SF_BOOT_GRAV_TIMEOUT_SEC')
 if timeout is not None and timeout > MAX_BOOT_TIMEOUT_SEC:
  bad.append(f'SF_BOOT_GRAV_TIMEOUT_SEC>{MAX_BOOT_TIMEOUT_SEC}')
 if timeout is not None:
  min_sec = p.get('SF_BOOT_GRAV_MIN_SEC', 0.0)
  hold_sec = p.get('SF_BOOT_GRAV_HOLD_SEC', 0.0)
  if timeout <= min_sec + hold_sec + MIN_MARGIN_SEC:
   bad.append('SF_BOOT_GRAV_TIMEOUT_SEC<=min+hold+margin')
 for key in p:
  if key in REJECTED_KEYS:
   bad.append(f'forbidden_key:{key}')
  if key.startswith(REJECTED_PREFIXES):
   bad.append(f'forbidden_key:{key}')
 return len(bad)==0, ';'.join(bad)

def parse(text):
 ds=RX['ds'].findall(text); ang=[m.groups() for m in RX['ang'].finditer(text)]; xyz=[m.groups() for m in RX['xyz'].finditer(text)]; r3d=[m.group(1) for m in RX['r3d'].finditer(text)]; accb=[m.groups() for m in RX['accb'].finditer(text)]; gates=[m.groups() for m in RX['gate'].finditer(text)]
 n=max(len(ds),len(ang),len(xyz),len(accb),1); out=[]
 for i in range(n):
  g=lambda arr,j,idx=float('nan'): float(arr[j][idx]) if j<len(arr) else float('nan')
  gp=(gates[i][0]=='1') if i<len(gates) else False; reason=(gates[i][1].strip() if i<len(gates) else 'missing_gate_status')
  out.append(dict(wave_dataset=ds[i] if i<len(ds) else f'dataset_{i}',roll_rms=g(ang,i,0),pitch_rms=g(ang,i,1),yaw_rms=g(ang,i,2),x_rms=g(xyz,i,0),y_rms=g(xyz,i,1),z_rms=g(xyz,i,2),rms_3d=float(r3d[i]) if i<len(r3d) else float('nan'),acc_bias_rms_3d=g(accb,i,3),quality_gate_pass=gp,fail_reason=reason))
 return out

def run_candidate(fam,cid,p,tier,seed,collect,stage):
 tdir=ROOT/'tests'/('kalman_ou_ii' if fam=='OU_II' else 'kalman_ou_iii'); bin_name='./kalman_ou_ii-sim' if fam=='OU_II' else './kalman_ou_iii-sim'
 env=os.environ.copy(); env['W3D_SEED']=str(seed); env['W3D_TIER']=tier
 if collect: env['W3D_COLLECT_ALL_GATES']='1'
 env.update({k:(f'{v:.6g}' if isinstance(v,float) else str(v)) for k,v in p.items()})
 pr=subprocess.run([bin_name],cwd=tdir,text=True,capture_output=True,env=env)
 rows=parse(pr.stdout+'\n'+pr.stderr)
 for r in rows:
  r.update({'family':fam,'candidate':cid,'seed':seed,'tier':tier,'stage':stage,'returncode':pr.returncode,'roll_pitch_rms_norm':math.hypot(r['roll_rms'],r['pitch_rms']),'xy_rms':math.hypot(r['x_rms'],r['y_rms'])})
  r.update(p)
  r['SF_BOOT_GRAV_TIMEOUT_SEC']=p.get('SF_BOOT_GRAV_TIMEOUT_SEC', float('nan'))
  r['timeout_valid']=int((not math.isnan(r['SF_BOOT_GRAV_TIMEOUT_SEC'])) and r['SF_BOOT_GRAV_TIMEOUT_SEC'] <= MAX_BOOT_TIMEOUT_SEC)
  r['rejected_before_run']=0
  r['reject_reason']=''
 return rows

def percentile(vals,p):
 if not vals:return float('inf')
 s=sorted(vals); i=(len(s)-1)*p; lo=int(i); hi=min(len(s)-1,lo+1); f=i-lo; return s[lo]*(1-f)+s[hi]*f

def score_rows(rows):
 by=defaultdict(list)
 for r in rows:
  if not r.get('rejected_before_run'): by[(r['family'],r['candidate'],r['seed'],r['tier'])].append(r)
 base={(r['family'],r['wave_dataset'],r['seed'],r['tier']):r for r in rows if r['candidate']=='baseline' and not r.get('rejected_before_run')}
 for _,it in by.items():
  rp=[]; z=[]; r3=[]; yaw=[]; acc=[]; any_fail=False
  for r in it:
   b=base.get((r['family'],r['wave_dataset'],r['seed'],r['tier']),r)
   rf=lambda x: x if math.isfinite(x) and abs(x)>1e-9 else 1e-9
   r['roll_pitch_ratio']=rf(r['roll_pitch_rms_norm'])/rf(b.get('roll_pitch_rms_norm', float('nan'))); r['z_ratio']=rf(r['z_rms'])/rf(b.get('z_rms', float('nan'))); r['rms3d_ratio']=rf(r['rms_3d'])/rf(b.get('rms_3d', float('nan'))); r['yaw_ratio']=rf(r['yaw_rms'])/rf(b.get('yaw_rms', float('nan'))); r['acc_bias_ratio']=rf(r['acc_bias_rms_3d'])/rf(b.get('acc_bias_rms_3d', float('nan')))
   rp.append(r['roll_pitch_ratio']); z.append(r['z_ratio']); r3.append(r['rms3d_ratio']); yaw.append(r['yaw_ratio']); acc.append(r['acc_bias_ratio']); any_fail = any_fail or (not r['quality_gate_pass'])
  s=4.0*percentile(rp,0.75)+2.0*max(rp)+2.0*max(0,max(z)-1.02)+2.0*max(0,max(r3)-1.02)+0.5*percentile(yaw,0.75)+0.5*percentile(acc,0.75)+(100.0 if any_fail else 0.0)
  for r in it: r['score']=s

def aggregate_candidates(rows, family, tier):
 out=[]
 by=defaultdict(list)
 for r in rows:
  if r.get('family')==family and r.get('tier')==tier:
   by[(r['candidate'])].append(r)
 for cand, it in by.items():
  valid = all((rr.get('returncode')==0 and rr.get('quality_gate_pass') and not rr.get('rejected_before_run')) for rr in it)
  params = {k:v for k,v in it[0].items() if re.match(r'^(SF_|OU_)', k)}
  obj={'family':family,'candidate':cand,'tier':tier,'params':params,'valid':valid,'n_rows':len(it),'stage':it[0].get('stage','unknown')}
  if valid:
   obj.update({
    'mean_score':sum(r['score'] for r in it)/len(it), 'max_score':max(r['score'] for r in it),
    'mean_rms3d':sum(r['rms_3d'] for r in it)/len(it), 'mean_z_rms':sum(r['z_rms'] for r in it)/len(it),
    'mean_roll_pitch_rms':sum(r['roll_pitch_rms_norm'] for r in it)/len(it), 'mean_yaw_rms':sum(r['yaw_rms'] for r in it)/len(it)
   })
  else:
   obj.update({'mean_score':float('inf'),'max_score':float('inf'),'mean_rms3d':float('inf'),'mean_z_rms':float('inf'),'mean_roll_pitch_rms':float('inf'),'mean_yaw_rms':float('inf')})
  out.append(obj)
 return out

def rank_key(a):
 return (a['mean_score'], a['max_score'], a['mean_rms3d'], a['mean_z_rms'], a['mean_roll_pitch_rms'], a['mean_yaw_rms'])

def eval_candidates(fam, cand_list, tier, collect, seed, stage):
 rows=[]
 for cid,p in cand_list:
  ok, reason = validate_candidate(p)
  if not ok:
   rows.append({'family':fam,'candidate':cid,'seed':seed,'tier':tier,'stage':stage,'wave_dataset':'rejected_before_run','quality_gate_pass':False,'fail_reason':'rejected_before_run','returncode':-1,'score':float('inf'),'rejected_before_run':1,'reject_reason':reason,**p})
   continue
  rows.extend(run_candidate(fam,cid,p,tier,seed,collect,stage))
 score_rows(rows)
 return rows

def sample_local_perturbations(base_params, ranges, rng, n):
 out=[]
 keys=[k for k in ranges if k in base_params]
 for i in range(n):
  p=dict(base_params)
  for k in keys:
   lo, hi = ranges[k]
   v = float(p[k])
   if lo > 0.0 and hi > 0.0:
    logv = math.log(v)
    spread = 0.18 * (math.log(hi)-math.log(lo))
    nv = math.exp(logv + rng.uniform(-spread, spread))
   else:
    spread = 0.18*(hi-lo)
    nv = v + rng.uniform(-spread, spread)
   p[k] = min(max(nv, lo), hi)
  if 'SF_BOOT_GRAV_MIN_SEC' in p and 'SF_BOOT_GRAV_HOLD_SEC' in p:
   apply_timeout_constraint(p, rng)
  out.append((f'local_{i:04d}', p))
 return out

def print_cpp_config(fam, result):
 params = result['params']
 print(f"\n=== BEST_CONFIG {fam} candidate={result['candidate']} tier={result['tier']} rows={result['n_rows']} ===")
 print(f"RMS_SUMMARY mean_3d={result['mean_rms3d']:.6g} mean_z={result['mean_z_rms']:.6g} mean_roll_pitch={result['mean_roll_pitch_rms']:.6g} mean_yaw={result['mean_yaw_rms']:.6g}")
 print("// C++ snippet:")
 print("struct TunedParams {")
 for k in sorted(params.keys()):
  v = params[k]
  if isinstance(v, float):
   print(f"  static constexpr float {k} = {v:.9g}f;")
  elif isinstance(v, int):
   print(f"  static constexpr int {k} = {v};")
 print("};")

def corr(xs, ys):
 n=len(xs)
 if n<2: return float('nan')
 mx=sum(xs)/n; my=sum(ys)/n
 vx=sum((x-mx)**2 for x in xs); vy=sum((y-my)**2 for y in ys)
 if vx<1e-14 or vy<1e-14: return float('nan')
 cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
 return cov / math.sqrt(vx*vy)

def print_param_sensitivity(rows, fam):
 valid=[r for r in rows if r.get('family')==fam and r.get('returncode')==0 and r.get('quality_gate_pass') and not r.get('rejected_before_run')]
 param_data=defaultdict(lambda: {'x':[], 'y':[]})
 for r in valid:
  for k,v in r.items():
   if re.match(r'^(SF_|OU_)', k) and isinstance(v,(float,int)) and math.isfinite(float(v)):
    param_data[k]['x'].append(float(v)); param_data[k]['y'].append(float(r.get('score', float('nan'))))
 print(f"\nPARAM_SENSITIVITY {fam}")
 print("param,n,min,max,corr_with_score")
 for k in sorted(param_data):
  xs=param_data[k]['x']; ys=param_data[k]['y']
  c=corr(xs,ys); mn=min(xs); mx=max(xs)
  print(f"{k},{len(xs)},{mn:.6g},{mx:.6g},{c if math.isfinite(c) else 'nan'}")
  if abs(mx-mn) < 1e-9:
   print(f"WARNING_PARAM_NO_VARIANCE {k}")

def write_csv(path, rows):
 keys=set()
 for r in rows: keys.update(r.keys())
 fields=sorted(keys)
 with open(path, 'w', newline='') as f:
  w=csv.DictWriter(f, fieldnames=fields)
  w.writeheader(); w.writerows(rows)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--mode',choices=['gravity','ou','full'],default='gravity'); ap.add_argument('--family',choices=['OU_II','OU_III','both'],default='both'); ap.add_argument('--samples',type=int,default=100); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--tier',choices=['quick','all','final'],default='quick'); ap.add_argument('--top-k',type=int,default=8); ap.add_argument('--collect-all-gates',action='store_true',default=False); ap.add_argument('--out-csv',default='')
 a=ap.parse_args(); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_ii'),'build'],check=True); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_iii'),'build'],check=True)
 fams=['OU_II','OU_III'] if a.family=='both' else [a.family]; all_rows=[]
 for fam in fams:
  rng=random.Random(a.seed+(0 if fam=='OU_II' else 1000000))
  # Stage A/B baseline+gravity
  grav_mc=[(f'grav_mc_{i:04d}',p) for i,p in enumerate(sample_params(GRAVITY_RANGES[fam],a.samples,rng,gravity=True),1)]
  stage_b=[('baseline',{})]+minimal_candidates()+explicit_practical_candidates()+grav_mc
  rows_b=eval_candidates(fam, stage_b, tier=a.tier, collect=(a.collect_all_gates or a.tier!='final'), seed=a.seed, stage='B')
  all_rows.extend(rows_b)
  agg_b=sorted([x for x in aggregate_candidates(rows_b, fam, a.tier) if x['valid']], key=rank_key)
  gravity_winners=agg_b[:max(1,a.top_k)]

  final_candidates=gravity_winners
  if a.mode in ('ou','full'):
   # Stage C OU on top of gravity winners
   ou_space=dict(OU_COMMON)
   if fam=='OU_II': ou_space.update(OU_II_ONLY)
   sampled_ou=sample_params(ou_space, a.samples, rng)
   stage_c=[]
   for gi,g in enumerate(gravity_winners,1):
    basep=g['params']
    for oi,ou in enumerate(sampled_ou,1):
      comb={**basep, **ou}
      stage_c.append((f'comb_g{gi:02d}_ou{oi:04d}', comb))
   rows_c=eval_candidates(fam, stage_c, tier=a.tier, collect=(a.collect_all_gates or a.tier!='final'), seed=a.seed, stage='C')
   all_rows.extend(rows_c)
   agg_c=sorted([x for x in aggregate_candidates(rows_c, fam, a.tier) if x['valid']], key=rank_key)
   top_c=agg_c[:max(1,a.top_k)]

   # Stage D local refinement
   local=[]
   local_ranges=dict(GRAVITY_RANGES[fam]); local_ranges.update(ou_space)
   for idx,c in enumerate(top_c,1):
    local.extend([(f'refine_t{idx:02d}_{cid}', p) for cid,p in sample_local_perturbations(c['params'], local_ranges, rng, max(1, a.samples//max(1,a.top_k)))])
   rows_d=eval_candidates(fam, local, tier=a.tier, collect=(a.collect_all_gates or a.tier!='final'), seed=a.seed, stage='D')
   all_rows.extend(rows_d)
   agg_d=sorted([x for x in aggregate_candidates(rows_d, fam, a.tier) if x['valid']], key=rank_key)
   final_candidates=(agg_d[:max(1,a.top_k)] or top_c or gravity_winners)

  # Stage E multi-seed final validation
  validate_tier = a.tier if a.tier != 'quick' else 'final'
  stage_e=[]
  for i,c in enumerate(final_candidates,1):
   for s in (a.seed, a.seed+1, a.seed+2):
    stage_e.append((f'final_{i:02d}_s{s}', c['params'], s))
  rows_e=[]
  for cid,p,s in stage_e:
   rows_e.extend(eval_candidates(fam, [(cid,p)], tier=validate_tier, collect=True, seed=s, stage='E'))
  all_rows.extend(rows_e)
  score_rows(all_rows)
  agg_e=sorted([x for x in aggregate_candidates(rows_e, fam, validate_tier) if x['valid']], key=rank_key)
  if not agg_e:
   print(f"\n=== BEST_CONFIG {fam} ===\nNO_VALID_CANDIDATES")
   continue
  best=agg_e[0]

  base=next((x for x in aggregate_candidates(rows_b, fam, a.tier) if x['candidate']=='baseline' and x['valid']), None)
  if base:
   ratio = best['mean_score']/base['mean_score'] if base['mean_score']>0 else float('inf')
   print(f"BASELINE_SCORE {base['mean_score']:.6g}")
   print(f"BEST_SCORE {best['mean_score']:.6g}")
   print(f"SCORE_RATIO {ratio:.6g}")
   if ratio > 0.99:
    print("NO_MEANINGFUL_IMPROVEMENT")
  best_rows=[r for r in rows_e if r['candidate']==best['candidate']]
  byds=defaultdict(list)
  for r in best_rows:
   byds[r['wave_dataset']].append(r.get('rms3d_ratio', float('nan')))
  print("PER_DATASET_RATIOS")
  for ds,v in sorted(byds.items()):
   print(f"{ds}: {sum(v)/len(v):.6g}")
  print(f"BEST_PARAMS {best['params']}")
  print_param_sensitivity(all_rows, fam)
  print_cpp_config(fam, best)

 if a.out_csv:
  write_csv(a.out_csv, all_rows)

if __name__=='__main__': main()
