#!/usr/bin/env python3
import os, re, csv, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CASES=[("baseline",{}),("racc_0_8",{"SF_RACC_WARMUP_STD":"0.8"}),("racc_1_2",{"SF_RACC_WARMUP_STD":"1.2"}),("racc_1_6",{"SF_RACC_WARMUP_STD":"1.6"}),("startup_A_300",{"SF_MAG_DELAY_SEC":"15","SF_MAG_GRAV_ALIGN_MAX_SIN":"0.07","SF_MAG_GRAV_ALIGN_HOLD_SEC":"2.0","SF_MAG_GRAV_ALIGN_LPF_TAU":"1.0","SF_MAG_TILT_FALLBACK_SEC":"30","SF_MAG_EXTREME_GYRO_DPS":"45","SF_BOOT_TILT_ACC_TAU":"2.5","SF_BOOT_GRAV_SLOW_TAU":"6.0","SF_BOOT_GRAV_ALIGN_MAX_SIN":"0.07","SF_BOOT_GRAV_HOLD_SEC":"2.0","SF_BOOT_GRAV_MIN_SEC":"8.0","SF_BOOT_GRAV_TIMEOUT_SEC":"15.0","SF_BOOT_GRAV_NORM_FRAC":"0.25","SF_ONLINE_TUNE_WARMUP_SEC":"10.0","SF_RACC_WARMUP_STD":"1.2","SF_MAG_MIN_SAMPLES":"300"}),("startup_A_600",{"SF_MAG_MIN_SAMPLES":"600","SF_MAG_DELAY_SEC":"15","SF_MAG_GRAV_ALIGN_MAX_SIN":"0.07","SF_MAG_GRAV_ALIGN_HOLD_SEC":"2.0","SF_MAG_GRAV_ALIGN_LPF_TAU":"1.0","SF_MAG_TILT_FALLBACK_SEC":"30","SF_MAG_EXTREME_GYRO_DPS":"45","SF_BOOT_TILT_ACC_TAU":"2.5","SF_BOOT_GRAV_SLOW_TAU":"6.0","SF_BOOT_GRAV_ALIGN_MAX_SIN":"0.07","SF_BOOT_GRAV_HOLD_SEC":"2.0","SF_BOOT_GRAV_MIN_SEC":"8.0","SF_BOOT_GRAV_TIMEOUT_SEC":"15.0","SF_BOOT_GRAV_NORM_FRAC":"0.25","SF_ONLINE_TUNE_WARMUP_SEC":"10.0","SF_RACC_WARMUP_STD":"1.2"})]
metric_re=re.compile(r"Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw=([0-9.eE+-]+)")
xyz_re=re.compile(r"XYZ RMS \(m\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+)")
def run_family(fam):
 tdir=ROOT/"tests"/("kalman_ou_ii" if fam=="OU_II" else "kalman_ou_iii")
 bin_name="./kalman_ou_ii-sim" if fam=="OU_II" else "./kalman_ou_iii-sim"; out=[]
 for cid,env_add in CASES:
  env=os.environ.copy(); env.update(env_add)
  p=subprocess.run([bin_name],cwd=tdir,env=env,text=True,capture_output=True)
  text=p.stdout+p.stderr; y=[float(m.group(3)) for m in metric_re.finditer(text)]; z=[float(m.group(3)) for m in xyz_re.finditer(text)]
  out.append(dict(family=fam,candidate=cid,returncode=p.returncode,yaw_rms_mean=(sum(y)/len(y) if y else 1e9),z_rms_mean=(sum(z)/len(z) if z else 1e9)))
 return out
subprocess.run(["make","-C",str(ROOT/"tests/kalman_ou_ii"),"build"],check=True)
subprocess.run(["make","-C",str(ROOT/"tests/kalman_ou_iii"),"build"],check=True)
rows=run_family("OU_II")+run_family("OU_III")
outp=ROOT/"tests"/"tuning_sweep_results.csv"
with outp.open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(outp)
