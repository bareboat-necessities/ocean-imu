"""Replay exported realized core operations; all numerical margins are diagnostic.

The high-precision oracle uses dense covariance algebra only on a short prefix.
The production word uses square-root covariance, three-row innovations and QR.
No fitted covariance, event selection or physical-history restart is permitted.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, solve_triangular, svd

from .factor_word import FactorWord, nuisance_residual
from .factor_certificates import singular_floor
from .matrix_certificates import ldlt


def process_blocks(q):
    groups=[list(range(6))]+[[6+i,9+i,12+i,15+i] for i in range(3)]+[list(range(18,21))]
    u=np.zeros_like(q); reconstructed=np.zeros_like(q)
    for idx in groups:
        block=q[np.ix_(idx,idx)]
        ldlt([[Fraction(float(x)) for x in row] for row in block])
        u[np.ix_(idx,idx)]=cholesky(block,lower=True)
        reconstructed[np.ix_(idx,idx)]=block
    if not np.array_equal(reconstructed,q):
        raise ValueError('cross-process blocks require a joint factor')
    return u


def high_precision_prefix(root,model,events,digits=70):
    """Independent short-prefix oracle, not a rigorous roundoff enclosure."""
    import mpmath as mp
    with mp.workdps(digits):
        def mat(a):
            return mp.matrix([[mp.mpf(float(v)) for v in row] for row in a])
        p=mat(root); c=mp.cholesky(p); t=mp.eye(len(root))
        f,q=mat(model['F']),mat(model['Q'])
        for event in events:
            if event['kind']=='prediction':
                p=f*p*f.T+q;t=f*t
            elif event['kind']=='correction':
                m=model[event['channel']];h,r=mat(m['H']),mat(m['R'])
                pc=p*h.T;s=h*pc+r;k=pc*(s**-1)
                a=mp.eye(len(root))-k*h
                p=a*p*a.T+k*r*k.T;t=a*t
        e=c.T*t.T*(p**-1)*t*c;e=(e+e.T)/2
        rho=mp.eigsy(e,eigvals_only=True)[-1,0]
        return {'digits':digits,'events':len(events),
                'rho_diagnostic':mp.nstr(rho,45),'delta_diagnostic':mp.nstr(1-rho,45),
                'rigorous_enclosure':False}


def replay(path,high_precision_events=60):
    rows=[json.loads(line) for line in Path(path).read_text().splitlines()]
    root=rows[0];model=rows[1];end=rows[-1];events=rows[2:-1]
    if root['kind']!='root' or model['kind']!='model' or end['kind']!='endpoint':
        raise ValueError('incomplete source export')
    p0=np.asarray(root['P'],float); f=np.asarray(model['F'],float)
    u=process_blocks(np.asarray(model['Q'],float));word=FactorWord(p0)
    checks={'gain_relative_max':0.,'innovation_relative_max':0.}
    e_hb=np.zeros((21,2));e_hb[2,0]=1;e_hb[5,1]=.02
    service=e_hb.copy();magrows=[];service_floors=[];sample=0
    short=FactorWord(p0)
    for index,event in enumerate(events):
        if event['kind']=='prediction':
            if sample and sample%200==0:
                values=svd(np.vstack(magrows),compute_uv=False)
                service_floors.append(float(values[-1]**2));magrows=[];service=e_hb.copy()
            sample+=1;word.prediction(f,u);service=f@service
            if index<high_precision_events:short.prediction(f,u)
        elif event['kind']=='correction':
            channel=event['channel'];m=model[channel];h=np.asarray(m['H'],float)
            actual_s=np.asarray(event['S_actual'],float);actual_k=np.asarray(event['gain_actual'],float)
            out=word.correction(h,m['R'])
            if index<high_precision_events:short.correction(h,m['R'])
            if channel=='mag':
                magrows.append(solve_triangular(cholesky(actual_s,lower=True),h@service,lower=True))
            service=service-actual_k@(h@service)
            for key,a,b in [('gain',out['gain'],actual_k),('innovation',out['innovation'],actual_s)]:
                error=float(np.linalg.norm(a-b,'fro')/max(np.linalg.norm(b,'fro'),1e-300))
                checks[key+'_relative_max']=max(checks[key+'_relative_max'],error)
        else:
            raise ValueError('unsupported exported event')
    if sample%200==0 and magrows:
        values=svd(np.vstack(magrows),compute_uv=False);service_floors.append(float(values[-1]**2))
    actual_p=np.asarray(end['P_actual'],float)
    checks['endpoint_covariance_relative']=float(np.linalg.norm(word.c@word.c.T-actual_p,'fro')/np.linalg.norm(actual_p,'fro'))
    scales=np.array([1]*3+[.02]*3+[2.4]*3+[18]*3+[132]*3+[4]*3+[1]*3)
    physical=solve_triangular(word.root_factor.T,word.loss.T,lower=False).T*scales
    nuisance=nuisance_residual(physical,[2,5],rank_tolerance=1e-12)
    ns=svd(nuisance['residual_factor'],compute_uv=False)
    hp=None
    if high_precision_events:
        hp=high_precision_prefix(p0,model,events[:high_precision_events])
        hp['binary64_rho_diagnostic']=short.snapshot()['rho_diagnostic']
        hp['absolute_parity_error']=abs(float(hp['rho_diagnostic'])-hp['binary64_rho_diagnostic'])
    x=solve_triangular(word.loss,np.eye(21),lower=False)
    matrix_bound=singular_floor(word.loss,x,factor_error_norm_upper=0)
    report={
        'qualification':'OU3_SHIPPING_CORE_FACTORIZED_WORD_DIAGNOSTIC_V1',
        'scope':root['scope'],'warmup_samples':root['warmup_samples'],
        'word_seconds':sample*root['dt'],'word':word.snapshot(),
        'shipping_double_parity':checks,'recorded_process_block_dimensions':[6,4,4,4,3],
        'recorded_process_blocks_spd_exact':True,
        'disjoint_service_windows':len(service_floors),
        'restricted_service_min_diagnostic':min(service_floors) if service_floors else None,
        'all_time_or_every_sliding_window_service_verified':False,
        'nuisance_rank_diagnostic':nuisance['nuisance_rank_diagnostic'],
        'heading_after_nuisance_min_diagnostic':float(ns[-1]**2),
        'computed_factor_matrix_only_bound':matrix_bound,
        'factor_matrix_bound_scope':'exact bound on the stored binary64 QR matrix only; whole-word factor error is not certified',
        'source_uniform_verified':False,'float32_transfer_verified':False,'theorem_closed':False,
    }
    if hp is not None:
        report['high_precision_prefix']=hp
    source=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
    report['shipping_header_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    return report


def main():
    parser=argparse.ArgumentParser();parser.add_argument('trace',type=Path)
    parser.add_argument('--output',type=Path);parser.add_argument('--high-precision-events',type=int,default=60)
    args=parser.parse_args()
    if not 0<=args.high_precision_events<=200:
        parser.error('high-precision prefix must contain 0 to 200 events')
    report=replay(args.trace,args.high_precision_events)
    rendered=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered,end='')
    # Diagnostic consistency gates, never a theorem promotion gate.
    return int(report['word']['identity_residual_fro']>1e-8 or
               max(report['shipping_double_parity'].values())>1e-7 or
               report.get('high_precision_prefix',{}).get('absolute_parity_error',0)>1e-11)


if __name__=='__main__':
    raise SystemExit(main())
