#!/usr/bin/env python3
import argparse, csv, math, os, random, re, subprocess, time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports' / 'results'
MAX_BOOT_TIMEOUT_SEC = 15.0
MIN_MARGIN_SEC = 1.0
MAX_MARGIN_SEC = 4.0
REJECTED_KEYS = {'OU_TAU_COEFF', 'OU_SIGMA_COEFF'}
REJECTED_PREFIXES = ('SF_MAG_',)

GRAVITY_RANGES = {
    'OU_II': {
        'SF_BOOT_TILT_ACC_TAU': (1.5, 2.8), 'SF_BOOT_GRAV_SLOW_TAU': (3.5, 6.0), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.12, 0.20),
        'SF_BOOT_GRAV_HOLD_SEC': (0.6, 1.5), 'SF_BOOT_GRAV_MIN_SEC': (4.0, 8.0), 'SF_BOOT_GRAV_NORM_FRAC': (0.30, 0.49),
        'SF_ONLINE_TUNE_WARMUP_SEC': (5.0, 10.0),
    },
    'OU_III': {
        'SF_BOOT_TILT_ACC_TAU': (1.5, 3.0), 'SF_BOOT_GRAV_SLOW_TAU': (3.5, 6.5), 'SF_BOOT_GRAV_ALIGN_MAX_SIN': (0.10, 0.20),
        'SF_BOOT_GRAV_HOLD_SEC': (0.6, 1.8), 'SF_BOOT_GRAV_MIN_SEC': (4.0, 8.5), 'SF_BOOT_GRAV_NORM_FRAC': (0.25, 0.49),
        'SF_ONLINE_TUNE_WARMUP_SEC': (5.0, 12.0),
    }
}
OU_COMMON = {
    'OU_ACC_NOISE_FLOOR_SIGMA': (0.10, 0.18),
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

def run_candidate(fam,cid,p,tier,seed,collect):
 tdir=ROOT/'tests'/('kalman_ou_ii' if fam=='OU_II' else 'kalman_ou_iii'); bin_name='./kalman_ou_ii-sim' if fam=='OU_II' else './kalman_ou_iii-sim'
 env=os.environ.copy(); env['W3D_SEED']=str(seed); env['W3D_TIER']=tier
 if collect: env['W3D_COLLECT_ALL_GATES']='1'
 env.update({k:(f'{v:.6g}' if isinstance(v,float) else str(v)) for k,v in p.items()})
 pr=subprocess.run([bin_name,'--nomag'],cwd=tdir,text=True,capture_output=True,env=env)
 rows=parse(pr.stdout+'\n'+pr.stderr)
 for r in rows:
  r.update({'family':fam,'candidate':cid,'seed':seed,'tier':tier,'mode':'gravity' if cid.startswith(('grav_','practical_','minimal_')) or cid in ('baseline',) else 'ou', 'returncode':pr.returncode,'roll_pitch_rms_norm':math.hypot(r['roll_rms'],r['pitch_rms']),'xy_rms':math.hypot(r['x_rms'],r['y_rms'])})
  r.update(p)
  r['SF_BOOT_GRAV_TIMEOUT_SEC']=p.get('SF_BOOT_GRAV_TIMEOUT_SEC', float('nan'))
  r['timeout_valid']=int((not math.isnan(r['SF_BOOT_GRAV_TIMEOUT_SEC'])) and r['SF_BOOT_GRAV_TIMEOUT_SEC'] <= MAX_BOOT_TIMEOUT_SEC)
  r['rejected_before_run']=0
  r['reject_reason']=''
 return rows

def progress_line(fam, mode, tier, i, total, cid, seed, rows_count, passed, failed, score, best_cand, best_score, start_ts, out_csv):
 elapsed = time.time() - start_ts
 eta = (elapsed / i) * (total - i) if i > 0 else 0
 fmt = lambda t: time.strftime('%H:%M:%S', time.gmtime(max(t, 0.0)))
 ts = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
 return f"{ts} [{fam}][{mode}][{tier}] {i}/{total} candidate={cid} seed={seed} rows={rows_count} pass={passed} fail={failed} score={score if score is not None else 'na'} best={best_cand} best_score={best_score if best_score is not None else 'na'} elapsed={fmt(elapsed)} eta={fmt(eta)} out={out_csv}"

def eval_set(fam, cand_list, tier='quick', collect=True, seed=42, out_csv=''):
 rows=[]; total=len(cand_list); start_ts=time.time(); prog_path=REPORTS/f'ou_sweep_progress_{fam.lower()}_{tier}_seed{seed}.log'
 best_score = None; best_cand = 'na'; pass_count=0; fail_count=0
 with prog_path.open('a') as prog:
  for i,(cid,p) in enumerate(cand_list,1):
   ok, reason = validate_candidate(p)
   if not ok:
    reject_row = {'family':fam,'candidate':cid,'seed':seed,'tier':tier,'mode':'gravity' if cid.startswith(('grav_','practical_','minimal_')) or cid in ('baseline',) else 'ou',
      'wave_dataset':'rejected_before_run','quality_gate_pass':False,'fail_reason':'rejected_before_run','returncode':-1,'score':float('inf'),'survived_tier':0,
      'SF_BOOT_GRAV_TIMEOUT_SEC':p.get('SF_BOOT_GRAV_TIMEOUT_SEC', float('nan')),'timeout_valid':0,'rejected_before_run':1,'reject_reason':reason}
    reject_row.update(p)
    rows.append(reject_row)
    fail_count += 1
   else:
    c_rows = run_candidate(fam,cid,p,tier,seed,collect)
    rows.extend(c_rows)
    cand_fail = any(not r['quality_gate_pass'] for r in c_rows)
    pass_count += 0 if cand_fail else 1
    fail_count += 1 if cand_fail else 0
   score_rows(rows)
   current = [r for r in rows if r['candidate']==cid and r['family']==fam and r['tier']==tier]
   cand_score = current[0].get('score') if current else None
   family_rows = [r for r in rows if r['family']==fam and r['tier']==tier and not r.get('rejected_before_run')]
   if family_rows:
    b = min(family_rows, key=lambda x: x.get('score', float('inf')))
    best_cand, best_score = b['candidate'], b.get('score')
   line = progress_line(fam, 'gravity' if cid.startswith(('grav_','practical_','minimal_')) or cid in ('baseline',) else 'ou', tier, i, total, cid, seed, len(current), pass_count, fail_count, cand_score, best_cand, best_score, start_ts, out_csv)
   print(line, flush=True)
   prog.write(line + '\n'); prog.flush()
 return rows

def percentile(vals,p):
 if not vals:return float('inf')
 s=sorted(vals); i=(len(s)-1)*p; lo=int(i); hi=min(len(s)-1,lo+1); f=i-lo; return s[lo]*(1-f)+s[hi]*f

def score_rows(rows):
 by=defaultdict(list)
 for r in rows:
  if not r.get('rejected_before_run'): by[(r['family'],r['candidate'],r['tier'])].append(r)
 base={(r['family'],r['wave_dataset'],r['seed'],r['tier']):r for r in rows if r['candidate']=='baseline' and not r.get('rejected_before_run')}
 for _,it in by.items():
  rp=[]; z=[]; r3=[]; yaw=[]; acc=[]; any_fail=False
  for r in it:
   b=base.get((r['family'],r['wave_dataset'],r['seed'],r['tier']),r)
   rf=lambda x: x if math.isfinite(x) and abs(x)>1e-9 else 1e-9
   r['baseline_roll_pitch_rms_norm']=b.get('roll_pitch_rms_norm', float('nan')); r['baseline_z_rms']=b.get('z_rms', float('nan')); r['baseline_rms_3d']=b.get('rms_3d', float('nan')); r['baseline_yaw_rms']=b.get('yaw_rms', float('nan'))
   r['roll_pitch_ratio']=rf(r['roll_pitch_rms_norm'])/rf(b.get('roll_pitch_rms_norm', float('nan'))); r['z_ratio']=rf(r['z_rms'])/rf(b.get('z_rms', float('nan'))); r['rms3d_ratio']=rf(r['rms_3d'])/rf(b.get('rms_3d', float('nan'))); r['yaw_ratio']=rf(r['yaw_rms'])/rf(b.get('yaw_rms', float('nan'))); r['acc_bias_ratio']=rf(r['acc_bias_rms_3d'])/rf(b.get('acc_bias_rms_3d', float('nan')))
   rp.append(r['roll_pitch_ratio']); z.append(r['z_ratio']); r3.append(r['rms3d_ratio']); yaw.append(r['yaw_ratio']); acc.append(r['acc_bias_ratio']); any_fail = any_fail or (not r['quality_gate_pass'])
  s=4.0*percentile(rp,0.75)+2.0*max(rp)+2.0*max(0,max(z)-1.02)+2.0*max(0,max(r3)-1.02)+0.5*percentile(yaw,0.75)+0.5*percentile(acc,0.75)+(100.0 if any_fail else 0.0)
  for r in it: r['score']=s; r['survived_tier']=int((not any_fail) and max(z)<=1.05 and max(r3)<=1.05)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--mode',choices=['gravity','ou','full'],default='gravity'); ap.add_argument('--family',choices=['OU_II','OU_III','both'],default='both'); ap.add_argument('--samples',type=int,default=100); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--tier',choices=['quick','all','final'],default='quick'); ap.add_argument('--top-k',type=int,default=8); ap.add_argument('--collect-all-gates',action='store_true',default=False)
 a=ap.parse_args(); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_ii'),'build'],check=True); subprocess.run(['make','-C',str(ROOT/'tests/kalman_ou_iii'),'build'],check=True)
 fams=['OU_II','OU_III'] if a.family=='both' else [a.family]; rows=[]
 REPORTS.mkdir(parents=True,exist_ok=True)
 out=REPORTS/f'ou_sweep_{a.mode}_{a.tier}_seed{a.seed}_n{a.samples}.csv'
 for fam in fams:
  rng=random.Random(a.seed+(0 if fam=='OU_II' else 1000000)); grav_mc=[(f'grav_mc_{i:04d}',p) for i,p in enumerate(sample_params(GRAVITY_RANGES[fam],a.samples,rng,gravity=True),1)]
  candidates=[('baseline',{})]+minimal_candidates()+explicit_practical_candidates()+grav_mc
  if a.mode in ('ou','full'):
   ou_space=dict(OU_COMMON)
   if fam=='OU_II': ou_space.update(OU_II_ONLY)
   candidates += [(f'ou_mc_{i:04d}',p) for i,p in enumerate(sample_params(ou_space,a.samples,rng),1)]
  rows.extend(eval_set(fam,candidates,tier=a.tier,collect=(a.collect_all_gates or a.tier!='final'),seed=a.seed,out_csv=str(out)))
 score_rows(rows)
 fields=sorted({k for r in rows for k in r.keys()})
 with out.open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
 print(out)

if __name__=='__main__': main()
