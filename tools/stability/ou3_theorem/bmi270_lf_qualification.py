#!/usr/bin/env python3
"""Analyze assembled AtomS3R/BMI270 low-frequency residual qualification captures."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np

G=9.80665
TX=60.0
THETA_MARGIN=0.0113349952420757

def design_fir(fs: float):
    # Deterministic symmetric windowed-sinc certification FIR.
    # Choose length for ~0.10-Hz transition; odd and capped by >=120 s support.
    n=int(math.ceil(4.0*fs/0.10))
    n=max(n,int(math.ceil(120*fs)))
    if n%2==0: n+=1
    m=np.arange(n)-(n-1)/2
    fc=.25
    h=2*fc/fs*np.sinc(2*fc/fs*m)
    h*=np.hamming(n)
    h/=h.sum()
    return h

def zero_phase_fir(x,h):
    # Symmetric FIR is zero-phase when evaluated centered offline; valid interior only.
    return np.vstack([np.convolve(x[:,j],h,mode="valid") for j in range(3)]).T

def read_csv(path):
    a=[]; g=[]; ar=[]; gr=[]; t=[]
    with open(path,newline="") as f:
        for r in csv.DictReader(f):
            t.append(float(r["t_s"]))
            a.append([float(r[k]) for k in ("ax_mps2","ay_mps2","az_mps2")])
            g.append([float(r[k]) for k in ("gx_rad_s","gy_rad_s","gz_rad_s")])
            ar.append([float(r[k]) for k in ("ax_ref_mps2","ay_ref_mps2","az_ref_mps2")])
            gr.append([float(r[k]) for k in ("gx_ref_rad_s","gy_ref_rad_s","gz_ref_rad_s")])
    return np.asarray(t),np.asarray(a)-np.asarray(ar),np.asarray(g)-np.asarray(gr)

def response(h,fs,f):
    n=np.arange(len(h))-(len(h)-1)/2
    return abs(np.sum(h*np.exp(-2j*np.pi*f*n/fs)))

def analyze(path,ua,ug):
    t,na,ng=read_csv(path)
    dt=np.diff(t)
    if len(dt)<2: raise ValueError("capture too short")
    fs=1/np.median(dt)
    if not (166.0 <= fs <= 251.0): raise ValueError(f"sample rate {fs:.3f} Hz outside qualified range")
    h=design_fir(fs)
    if len(t)<len(h)+10: raise ValueError(f"capture needs > {len(h)/fs:.1f} s for L_X interior")
    la=zero_phase_fir(na,h); lg=zero_phase_fir(ng,h)
    ea=float(np.linalg.norm(la,axis=1).max()+ua)
    eg=float(np.linalg.norm(lg,axis=1).max()+ug)
    charge=2*math.asin(min(1.0,ea/G))+TX*eg
    sha=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return {"file":str(path),"sha256":sha,"fs_hz":fs,"fir_taps":len(h),
            "gain_0p15":response(h,fs,.15),"gain_0p30":response(h,fs,.30),
            "eps_a_LF_mps2":ea,"eps_g_LF_rad_s":eg,
            "angular_charge_rad":charge,"margin_rad":THETA_MARGIN-charge,
            "pass":charge < THETA_MARGIN}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("captures",nargs="+")
    ap.add_argument("--u-a-ref",type=float,required=True)
    ap.add_argument("--u-g-ref",type=float,required=True)
    ap.add_argument("--output",required=True)
    ns=ap.parse_args()
    rows=[analyze(p,ns.u_a_ref,ns.u_g_ref) for p in ns.captures]
    ea=max(r["eps_a_LF_mps2"] for r in rows); eg=max(r["eps_g_LF_rad_s"] for r in rows)
    charge=2*math.asin(min(1.0,ea/G))+TX*eg
    out={"kind":"ASSEMBLED_BMI270_LF_RESIDUAL_QUALIFICATION",
         "qualified":charge<THETA_MARGIN,
         "T_X_s":TX,"theta_X_deg":1.0,
         "eps_a_LF_mps2":ea,"eps_g_LF_rad_s":eg,
         "joint_angular_charge_rad":charge,
         "allowed_charge_rad":THETA_MARGIN,
         "margin_rad":THETA_MARGIN-charge,
         "reference_uncertainty":{"u_a_mps2":ns.u_a_ref,"u_g_rad_s":ns.u_g_ref},
         "captures":rows}
    Path(ns.output).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
    raise SystemExit(0 if out["qualified"] else 2)
if __name__=="__main__": main()
