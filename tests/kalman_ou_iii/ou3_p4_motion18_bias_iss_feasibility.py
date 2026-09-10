#!/usr/bin/env python3
"""Point small-gain feasibility for the 18-state motion / bounded-bias P4 split.

This is deliberately non-promoting.  It consumes the same complete 21-state
shipping word/covariance trace as the motion18 backbone diagnostic and partitions
one composed word only after all shipping events have been multiplied:

    [x_m^+]   [ A  B ] [x_m]
    [e_b^+] = [ C  D ] [e_b]

where x_m is the 18-state motion coordinate and e_b the accelerometer-bias
coordinate.  The declared theorem does not require e_b to contract; projection
and the physical BIAS0/1/2 recurrence provide its compact/supply bound.  What is
needed is a strict homogeneous motion backbone A plus a finite gain from bounded
bias into motion.

With moving covariance metric V_k=x_m' P_k^-1 x_m, define

    rho0 = || P1^-1/2 A P0^1/2 ||_2^2,
    g2   = || P1^-1/2 B ||_2^2.

For rho=(1+rho0)/2 and eta=rho/rho0-1, Young's inequality gives

    V1 <= rho V0 + gamma_b ||e_b||^2,
    gamma_b = (1+1/eta) g2.

This is only a point feasibility test of the master inequality.  Production
proof must replace point matrices by same-history source-uniform outward cells,
charge BIAS-family driver/tau-mismatch supplies at every literal prefix, include
nonlinear/reset/projection graphs and binary32 ISS, and prove a retained tube.
"""
from __future__ import annotations

import argparse,json,math
from pathlib import Path
import numpy as np

import ou3_p4_complete_brmm_word_feasibility as WORD


def _sym(A): return 0.5*(A+A.T)

def _sqrt_psd(P):
    w,V=np.linalg.eigh(_sym(P))
    if float(np.min(w))<=0: raise RuntimeError('metric covariance block is not SPD')
    return (V*np.sqrt(w))@V.T

def _inv_sqrt_psd(P):
    w,V=np.linalg.eigh(_sym(P))
    if float(np.min(w))<=0: raise RuntimeError('metric covariance block is not SPD')
    return (V*(1.0/np.sqrt(w)))@V.T

def _gain(A,B,P0,P1):
    S0=_sqrt_psd(P0); J1=_inv_sqrt_psd(P1)
    T=J1@A@S0; G=J1@B
    rho0=float(np.linalg.norm(T,2)**2)
    g2=float(np.linalg.norm(G,2)**2)
    if not (math.isfinite(rho0) and math.isfinite(g2) and rho0>=0 and g2>=0):
        raise RuntimeError('non-finite motion/bias gain')
    if not rho0<1.0:
        return {'rho0':rho0,'strict_backbone':False,'g2':g2,'rho_target':None,'eta':None,'gamma_bias':None}
    rho=0.5*(1.0+rho0)
    if rho0==0.0:
        # No cross term from homogeneous motion; direct bias term alone remains.
        return {'rho0':rho0,'strict_backbone':True,'g2':g2,'rho_target':rho,'eta':math.inf,'gamma_bias':g2}
    eta=rho/rho0-1.0
    gamma=(1.0+1.0/eta)*g2
    return {'rho0':rho0,'strict_backbone':True,'g2':g2,'rho_target':rho,'eta':eta,'gamma_bias':gamma}

def analyze(map_path:Path,cov_path:Path,motion_json:Path):
    ms,maps=WORD.read_maps(map_path);cs,covs=WORD.read_covariances(cov_path)
    motion=json.loads(motion_json.read_text())
    if ms!=cs:raise RuntimeError('map/cov stride mismatch')
    out={}
    for mode in ('H18','A21'):
        worst=motion['modes'][mode]['worst']
        if worst is None:raise RuntimeError(mode+' has no motion word')
        i=int(worst['record']);m=maps[i];c=covs[i]
        A=m['M'][:18,:18];B=m['M'][:18,18:21]
        g=_gain(A,B,c['P0'][:18,:18],c['P1'][:18,:18])
        # Cross-check that we analyzed the exact same controlling motion block.
        if abs(g['rho0']-float(worst['rho_linear']))>2e-10*max(1.0,abs(float(worst['rho_linear']))):
            raise RuntimeError(f'{mode} rho mismatch {g["rho0"]} != {worst["rho_linear"]}')
        bnorm=float(np.linalg.norm(B,2))
        g.update({'record':i,'t0':float(m['t0']),'t1':float(m['t1']),
                  'raw_B_spectral_norm':bnorm,
                  'finite_bias_gain':g['gamma_bias'] is not None and math.isfinite(float(g['gamma_bias'])),
                  'bias_root_is_supply_not_homogeneous_motion_state':True})
        out[mode]=g
    return {
      'qualification':'NON_PROMOTING_OU3_P4_MOTION18_BOUNDED_BIAS_SMALL_GAIN_V1',
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','point_same_history_diagnostic_only':True,
      'P4_promoted':False,'filter_changed':False,'declared_domain_changed':False,
      'complete_21_state_word_partitioned_only_after_composition':True,
      'moving_covariance_motion_metric_used':True,
      'young_inequality_target':'rho=(1+rho0)/2; gamma_b=(1+1/eta)*||P1^-1/2 B||^2',
      'BIAS0_BIAS1_BIAS2_source_uniform_supply_not_closed_here':True,
      'nonlinear_finite_amplitude_not_closed_here':True,
      'every_prefix_retention_not_closed_here':True,
      'modes':out,
      'next_obligation':'if nonlinear motion18 shadow stays strict, replace this point A/B partition by source-uniform complete-word/prefix cells and absorb the existing all-family bias supply plus BRMM and binary32 channels'}
def validate(d):
    f=[]
    for k in ('point_same_history_diagnostic_only','complete_21_state_word_partitioned_only_after_composition','moving_covariance_motion_metric_used','BIAS0_BIAS1_BIAS2_source_uniform_supply_not_closed_here','nonlinear_finite_amplitude_not_closed_here','every_prefix_retention_not_closed_here'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('P4_promoted','filter_changed','declared_domain_changed'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode in ('H18','A21'):
        m=d.get('modes',{}).get(mode,{})
        if m.get('strict_backbone') is not True:f.append(mode+' backbone not strict')
        if m.get('finite_bias_gain') is not True:f.append(mode+' bias gain not finite')
        for k in ('rho0','rho_target','g2','gamma_bias'):
            x=m.get(k)
            if x is None or not math.isfinite(float(x)):f.append(mode+'.'+k+' invalid')
        if m.get('rho_target') is not None and not(float(m['rho0'])<float(m['rho_target'])<1.0):f.append(mode+' rho target not between rho0 and one')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--map',type=Path,required=True);ap.add_argument('--cov',type=Path,required=True);ap.add_argument('--motion-json',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=analyze(a.map,a.cov,a.motion_json);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'H18':d['modes']['H18'],'A21':d['modes']['A21'],'failures':f},indent=2,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
