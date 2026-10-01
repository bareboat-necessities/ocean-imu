#!/usr/bin/env python3
"""Finite-capture LF diagnostic (legacy filename), NEVER IMU qualification.

The former hard-coded 1-degree/60-second/FIR gate was not a proved general
physical gauge inequality. It is retired. The current contract is in
imu_bias.py/imu_temporal.py: supplied slow/fast histories and both all-window
accumulation profiles. LF magnitudes of a finite total-error record cannot
supply that decomposition or certify an all-time continuation.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np


def design_fir(fs: float, cutoff_hz: float, transition_hz: float):
    """Explicit user-supplied diagnostic filter, with no proof-selected cutoff."""
    if (not all(math.isfinite(x) and x>0 for x in (fs,cutoff_hz,transition_hz))
            or cutoff_hz+transition_hz/2>=fs/2):
        raise ValueError('positive finite filter settings below Nyquist required')
    n=max(3,int(math.ceil(4*fs/transition_hz)))
    if n%2==0:n+=1
    m=np.arange(n)-(n-1)/2
    h=2*cutoff_hz/fs*np.sinc(2*cutoff_hz/fs*m)*np.hamming(n)
    return h/h.sum()


def zero_phase_fir(x,h):
    return np.column_stack([np.convolve(x[:,j],h,mode='valid') for j in range(3)])


def read_csv(path):
    arrays=[[] for _ in range(5)]
    with open(path,newline='') as f:
        for r in csv.DictReader(f):
            arrays[0].append(float(r['t_s']))
            for out,keys in zip(arrays[1:],(
                ('ax_mps2','ay_mps2','az_mps2'),('gx_rad_s','gy_rad_s','gz_rad_s'),
                ('ax_ref_mps2','ay_ref_mps2','az_ref_mps2'),('gx_ref_rad_s','gy_ref_rad_s','gz_ref_rad_s'))):
                out.append([float(r[k]) for k in keys])
    t,a,g,ar,gr=map(np.asarray,arrays)
    if len(t)<3 or not all(np.isfinite(x).all() for x in (t,a,g,ar,gr)):
        raise ValueError('finite referenced samples required')
    return t,a-ar,g-gr


def response(h,fs,f):
    n=np.arange(len(h))-(len(h)-1)/2
    return abs(np.sum(h*np.exp(-2j*np.pi*f*n/fs)))


def analyze(path,ua,ug,*,cutoff_hz,transition_hz):
    if not all(math.isfinite(x) and x>=0 for x in (ua,ug)):
        raise ValueError('finite nonnegative reference uncertainty required')
    t,na,ng=read_csv(path);dt=np.diff(t)
    if np.any(dt<=0):raise ValueError('strictly increasing delivered timestamps required')
    # This legacy discrete FIR diagnostic has no nonuniform-time transfer proof.
    if not np.allclose(dt,np.median(dt),rtol=1e-5,atol=1e-10):
        raise ValueError('uniform timestamps required by this FIR diagnostic; do not silently resample')
    fs=1/np.median(dt);h=design_fir(fs,cutoff_hz,transition_hz)
    if len(t)<len(h):raise ValueError('capture shorter than diagnostic filter support')
    la,lg=zero_phase_fir(na,h),zero_phase_fir(ng,h)
    # For a signed FIR, deterministic reference error pays its L1 gain.
    gain=float(np.abs(h).sum())
    return {'kind':'FINITE_TOTAL_ERROR_LF_DIAGNOSTIC_ONLY','qualified':False,
        'file':str(path),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        'fs_hz':float(fs),'fir_taps':len(h),'cutoff_hz':cutoff_hz,'transition_hz':transition_hz,
        'reference_error_L1_gain':gain,
        'eps_a_LF_mps2':float(np.linalg.norm(la,axis=1).max()+gain*ua),
        'eps_g_LF_rad_s':float(np.linalg.norm(lg,axis=1).max()+gain*ug),
        'slow_fast_decomposition_certified':False,'fast_temporal_qualification':'OPEN',
        'all_time_continuation_certified':False,'motion_threshold_inferred':False}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('captures',nargs='+');ap.add_argument('--u-a-ref',type=float,required=True)
    ap.add_argument('--u-g-ref',type=float,required=True)
    ap.add_argument('--cutoff-hz',type=float,required=True)
    ap.add_argument('--transition-hz',type=float,required=True)
    ap.add_argument('--output',required=True);ns=ap.parse_args()
    rows=[analyze(p,ns.u_a_ref,ns.u_g_ref,cutoff_hz=ns.cutoff_hz,transition_hz=ns.transition_hz)
          for p in ns.captures]
    out={'kind':'FINITE_LF_DIAGNOSTICS_ONLY','qualified':False,'captures':rows,
         'controlling_IMU_model':'SLOW+FAST; numerical temporal qualification OPEN'}
    Path(ns.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0  # Successful diagnostic, NOT successful physical qualification.


if __name__=='__main__':raise SystemExit(main())
