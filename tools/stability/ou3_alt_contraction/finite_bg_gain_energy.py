"""Exact same-history gyro-bias gain-energy accounting before any rho claim.

The shipping covariance reset and a_w floor leave the b_g marginal unchanged.
For each actual gain/innovation event, charge the positive defect in

    tr(K_bg R K_bg') <= tr(P_bg_before - P_bg_after).

This does NOT assume that the floating-point Joseph update is PSD-decreasing.
The defect is computed from the same stored operands; covariance prediction
increments and mean injection roundoff are retained. Telescoping and weighted
Cauchy--Schwarz then give a finite-prefix bias bound. A native trace supplies
only a per-history diagnostic, never uniform defect or innovation bounds.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import os
import subprocess

from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import shipping_finite_identity as NATIVE

QUALIFICATION='OU3_ALT_SAME_HISTORY_BG_GAIN_ENERGY_V1'
ROOT=Path(__file__).resolve().parents[3]
SOURCE=ROOT/'tests/ou3_alt_contraction/live_covariance_feasibility.cpp'
INNOVATIONS=(12,22,32)


def trace(p): return sum(p[i][i] for i in range(3))
def matrix(xs): return tuple(tuple(xs[3*i+j] for j in range(3)) for i in range(3))
def norm1(xs): return sum(map(abs,xs))


@dataclass(frozen=True)
class Row:
    step:int
    kind:int
    covariance:tuple
    bias:tuple
    gain:tuple|None=None
    noise:tuple|None=None
    residual:tuple|None=None


def rows(path):
    out=[]
    for line in Path(path).read_text().splitlines():
        values=tuple(map(int,line.split()))
        if len(values)<2: raise ValueError('incomplete bias event row')
        step,kind=values[:2]
        if len(values)!=(35 if kind in INNOVATIONS else 14):
            raise ValueError('bias event row schema mismatch')
        v=tuple(M.value(x) for x in values[2:])
        out.append(Row(step,kind,matrix(v[:9]),v[9:12],
                       matrix(v[12:21]) if kind in INNOVATIONS else None,
                       matrix(v[21:30]) if kind in INNOVATIONS else None,
                       v[30:33] if kind in INNOVATIONS else None))
    return tuple(out)


def account(events,*,steps=600):
    """Algebraic bound on one correlated trace, with no PSD replacement.

    Every predecessor is checked, including prediction, post-mean, post-Joseph
    and sample-boundary records. Diagonal positive R is the literal selected
    dormant-guard/native profile. General correlated R requires an exact SPD
    inverse witness and is deliberately not supplied by this diagnostic.
    """
    events=tuple(events)
    if not events or events[0].kind!=0 or events[0].step!=0:
        raise ValueError('actual startup/release seed row required')
    p,b=events[0].covariance,events[0].bias
    p0,b0=p,b
    energy=innovation_energy=defect=prediction=mean_roundoff=F(0)
    gain_events=completed=0
    largest_defect=F(0)
    max_bias_norm2=sum(x*x for x in b)
    cursor=1
    for step in range(1,steps+1):
        if cursor+1>=len(events): raise ValueError('incomplete prediction')
        before,after=events[cursor:cursor+2];cursor+=2
        if (before.step,before.kind,after.step,after.kind)!=(step,1,step,2):
            raise ValueError('actual prediction pair missing')
        if (before.covariance,before.bias)!=(p,b) or after.bias!=b:
            raise ValueError('prediction detached from carried bias/covariance')
        prediction+=trace(after.covariance)-trace(p)
        p=after.covariance
        while cursor<len(events) and events[cursor].kind!=100:
            if cursor+2>=len(events): raise ValueError('incomplete Joseph triple')
            inn,mean,cov=events[cursor:cursor+3];cursor+=3
            if inn.kind not in INNOVATIONS or (mean.kind,cov.kind)!=(inn.kind+1,inn.kind+2):
                raise ValueError('gain/mean/Joseph ancestry missing')
            if any(x.step!=step for x in (inn,mean,cov)):
                raise ValueError('gain event moved to different IMU prefix')
            if (inn.covariance,inn.bias)!=(p,b) or mean.covariance!=p or cov.bias!=mean.bias:
                raise ValueError('gain event detached from carried state')
            K,R,r=inn.gain,inn.noise,inn.residual
            if K is None or R is None or r is None: raise ValueError('same-event gain/noise/innovation required')
            if any(R[i][j] for i in range(3) for j in range(3) if i!=j) or any(R[i][i]<=0 for i in range(3)):
                raise ValueError('literal positive diagonal measurement covariance required')
            e=sum(K[i][j]**2*R[j][j] for i in range(3) for j in range(3))
            v=sum(r[j]**2/R[j][j] for j in range(3))
            d=max(F(0),e-(trace(p)-trace(cov.covariance)))
            kb=tuple(sum(K[i][j]*r[j] for j in range(3)) for i in range(3))
            mean_roundoff+=norm1(tuple(mean.bias[i]-b[i]-kb[i] for i in range(3)))
            energy+=e;innovation_energy+=v;defect+=d;largest_defect=max(largest_defect,d)
            p,b=cov.covariance,cov.bias;gain_events+=1
            max_bias_norm2=max(max_bias_norm2,sum(x*x for x in b))
        if cursor>=len(events): raise ValueError('complete IMU boundary required')
        end=events[cursor];cursor+=1
        if (end.step,end.kind)!=(step,100) or (end.covariance,end.bias)!=(p,b):
            raise ValueError('unaccounted bias/covariance change in reset/floor/suffix')
        completed+=1
    if cursor!=len(events): raise ValueError('trace has an unaccounted tail')
    exact_budget=trace(p0)+prediction+defect-trace(p)
    if not energy<=exact_budget: raise AssertionError('gain-energy telescoping identity failed')
    # Retain final trace explicitly; no PSD theorem is inferred from a trace.
    # sqrt(E*V) is represented by its exact squared upper. L1 is a safe norm
    # upper for initial bias and accumulated machine injection defects.
    additive=norm1(b0)+mean_roundoff
    cauchy_square=exact_budget*innovation_energy
    if min(exact_budget,cauchy_square)<0: raise AssertionError('negative energy budget')
    measured_norm2=sum(x*x for x in b)
    # The energy budget is nondecreasing: each Joseph contributes drop+D>=e,
    # prediction cancels its own exact trace increment. Both innovation energy
    # and the accumulated L1 injection charge are nondecreasing too. The final
    # bound therefore covers every recorded prefix, not just its final bias.
    if max_bias_norm2>additive*additive+cauchy_square:
        cross=max_bias_norm2-additive*additive-cauchy_square
        if cross*cross>4*additive*additive*cauchy_square:
            raise AssertionError('weighted Cauchy mean bound failed')
    return {'qualification':QUALIFICATION,'steps':completed,'gain_events':gain_events,
            'initial_bg_trace':trace(p0),'final_bg_trace':trace(p),
            'actual_prediction_trace_increment':prediction,
            'gain_energy':energy,'innovation_energy':innovation_energy,
            'Joseph_positive_defect_sum':defect,'largest_Joseph_positive_defect':largest_defect,
            'telescoping_energy_budget':exact_budget,'bias_norm_additive_upper':additive,
            'bias_norm_radical_square_upper':cauchy_square,
            'final_bias_norm_squared':measured_norm2,'max_recorded_bias_norm_squared':max_bias_norm2,
            'every_recorded_bias_prefix_covered_by_same_energy_bound':True,
            'same_history_telescoping_and_Cauchy_bound_closed':True,
            'floating_Joseph_PSD_nonincrease_assumed':False,
            'sample_maxima_are_source_uniform_ceilings':False,
            'source_uniform_defect_and_innovation_bounds_closed':False,
            'all_event_arithmetic_source_uniform':False,'storage_search_allowed':False}


def run_native(work,*,mode='H',scenario=0):
    """Run a source-generated 600-event native feasibility word.

    The source schedules are bounded residuals beside the stationary physical
    reference. The three selected schedules keep the actual guard dormant;
    this check is post-run regression evidence, not an admission definition.
    """
    if mode not in ('H','A') or scenario not in (0,1,2): raise ValueError('unknown diagnostic schedule')
    work=Path(work);work.mkdir(parents=True,exist_ok=True)
    overlay=work/'include';manifest=NATIVE.make_overlay(overlay)
    eigen=Path(os.environ.get('OU3_ALT_TARGET_EIGEN_INCLUDE_DIR','/tmp/ou3-target/eigen/Eigen-0.3.2/ArduinoEigen'))
    if not (eigen/'Eigen/Dense').is_file(): eigen=NATIVE.eigen_include()
    exe=work/'probe'
    fingerprint=hashlib.sha256((SOURCE.read_text()+json.dumps(manifest,sort_keys=True)+str(eigen)).encode()).hexdigest()
    receipt=work/'probe.sha256'
    if not exe.exists() or not receipt.exists() or receipt.read_text()!=fingerprint:
        command=[os.environ.get('CXX','g++'),'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                 '-DEIGEN_DONT_VECTORIZE','-I'+str(overlay),'-I'+str(eigen),'-I'+str(ROOT/'src'),
                 str(SOURCE),'-o',str(exe)]
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=180)
        receipt.write_text(fingerprint)
    path=work/f'{mode}-{scenario}.txt'
    proc=subprocess.run([str(exe),mode,str(scenario),str(path)],check=True,capture_output=True,text=True,timeout=60)
    lines=[line.split() for line in proc.stdout.splitlines()]
    if len(lines)!=2 or lines[0][0]!='seed' or lines[1][0]!='run': raise ValueError('native report schema drift')
    seed,run=lines
    if run[4]!='1' or F(run[5])!=0: raise ArithmeticError('native word nonfinite or guard not dormant')
    names=('P_spectral_max','P_eigenvalue_min','bg_norm_max','aw_norm_max','gain_Frobenius_max',
           'innovation_norm_max','S_eigenvalue_min','injection_angle_norm_max','reset_spectral_ratio_max',
           'reset_ratio_step','reset_ratio_injection_angle','Joseph_Loewner_increase_max','Pbg_spectral_max')
    measurements=dict(zip(names,map(float,run[6:])))
    if len(run[6:])!=len(names): raise ValueError('native metric schema drift')
    report=account(rows(path))
    report['native']={'mode':mode,'scenario':scenario,'startup_ticks':int(seed[3]),'live_tick':int(seed[4]),
                      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                      'source_manifest':manifest,'report_lines':proc.stdout.splitlines(),'measurements':measurements,
                      'raw_packet_caps_checked_before_updates':True,
                      'guard_dormant_in_this_trace':True,'rho_or_universal_cover_inferred':False}
    return report


def summary(report):
    fields=('initial_bg_trace','final_bg_trace','actual_prediction_trace_increment',
            'gain_energy','innovation_energy','Joseph_positive_defect_sum',
            'largest_Joseph_positive_defect','telescoping_energy_budget',
            'bias_norm_additive_upper','final_bias_norm_squared')
    return {'qualification':QUALIFICATION, 'mode':report['native']['mode'],
            'scenario':report['native']['scenario'],'steps':report['steps'],
            'gain_events':report['gain_events'],
            **{key:float(report[key]) for key in fields},
            'bias_norm_bound_approx':float(report['bias_norm_additive_upper'])+
                math.sqrt(float(report['bias_norm_radical_square_upper'])),
            'native':report['native']['measurements'],
            'source_uniform_defect_and_innovation_bounds_closed':False,
            'all_event_arithmetic_source_uniform':False,'storage_search_allowed':False}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--mode',choices=('H','A'),default='H')
    parser.add_argument('--scenario',type=int,choices=(0,1,2),default=0)
    args=parser.parse_args()
    print(json.dumps(summary(run_native(args.work,mode=args.mode,scenario=args.scenario)),indent=2))
