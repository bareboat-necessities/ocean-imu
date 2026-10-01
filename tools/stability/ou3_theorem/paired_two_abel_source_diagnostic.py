"""Carried paired physical-source diagnostic; never source-uniform promotion.

Builds the causal terminal-AW adjoint on the literal exported word, retains the
aligned physical p/v/S/a truth, and evaluates the SAME-history physical source
response by paired finite differences.  This is a feasibility diagnostic for
G_phys,W; the theorem still needs an interval/factor enclosure.
"""
from pathlib import Path
import argparse, json, math, subprocess, tempfile
import numpy as np
from .ag_readout_source_diagnostic import REPO, HEADER, instrument, validate_physical_lift
from .aw_port_17s_diagnostic import normalized_events

DRIVER=REPO/'tools/stability/ag_readout_source.cpp'
RESERVE=.003056118

def export(eigen, heading):
    source=(REPO/HEADER).read_text()
    with tempfile.TemporaryDirectory(prefix='ou3-paired-physical-') as d:
        d=Path(d); inc=d/'kalman_ou_iii'; inc.mkdir()
        (inc/HEADER.name).write_text(instrument(source))
        exe=d/'probe'
        subprocess.run(['g++','-O2','-std=c++20','-I'+str(d),'-I'+str(REPO/'src'),
                        '-isystem',str(eigen),str(DRIVER),'-o',str(exe)],check=True)
        return json.loads(subprocess.check_output([str(exe),heading],text=True))

def backward_aw_rows(events, axis):
    y=np.zeros(21); y[15+axis]=1.0
    rows=[None]*(len(events)+1); rows[-1]=y.copy()
    for i in range(len(events)-1,-1,-1):
        e=events[i]
        if e['kind']=='correction': y=y@(np.eye(21)-e['K']@e['H'])
        elif e['kind']=='prediction': y=y@e['F']
        else: y=y@e['G']
        rows[i]=y.copy()
    return rows

def physical_series(trace):
    # One physical sample per normalized event. sync_completion inherits the
    # same literal timestamp/primitive as the completed prediction boundary.
    raw=[]
    pending=None
    for e in trace['events']:
        if e['kind']=='prediction':
            pending=e
            raw.append(e)
        elif e['kind']=='sync':
            continue
        elif e['kind']=='sync_completion':
            raw.append(e)
            pending=None
        else:
            raw.append(e)
    if pending is not None: raise ValueError('unfinished physical prediction')
    return raw

def analyze(trace):
    validate_physical_lift(trace)
    events=normalized_events(trace)
    phys=physical_series(trace)
    if len(events)!=len(phys): raise ValueError('literal/physical event alignment mismatch')
    # Exact carried output due to the physical truth is reconstructed by
    # operation-local affine source increments: prediction uses the change of
    # the physical LIN primitive under the exported F_LL; S and accelerometer
    # corrections use their literal K times the physical innovation component.
    # Magnetic physical terms are deliberately excluded here: this diagnostic
    # targets the p/v/S/a paired bridge only.
    total=np.zeros(3)
    per_axis=[]
    for axis in range(3):
        rows=backward_aw_rows(events,axis)
        val=0.0; proc=0.0; corr=0.0
        prev=None
        for i,(e,pe) in enumerate(zip(events,phys)):
            x=np.concatenate([
                np.array(pe['physical_v'],float).reshape(3),
                np.array(pe['physical_p'],float).reshape(3),
                np.array(pe['physical_S'],float).reshape(3),
                np.array(pe['physical_a'],float).reshape(3)])
            if prev is None: prev=x.copy()
            if e['kind']=='prediction':
                f=e['F'][6:18,6:18]
                d=x-f@prev
                term=float(rows[i+1][6:18]@d); proc+=term; val+=term
                prev=x.copy()
            elif e['kind']=='correction' and pe.get('sensor') in ('acc','S'):
                k=e['K']
                if pe['sensor']=='S':
                    rphys=-x[6:9]
                else:
                    # Translational physical specific-force contribution in
                    # the exported world fixture; attitude/gravity terms belong
                    # to the separate nonlinear/joint-rotation bridge.
                    rphys=np.array(pe['physical_a'],float).reshape(3)
                term=float(rows[i+1]@k@rphys); corr+=term; val+=term
        total[axis]=val
        per_axis.append({'axis':axis,'paired_process':proc,'paired_acc_S':corr,
                         'paired_sum':val,'unpaired_abs_sum':abs(proc)+abs(corr)})
    norm=float(np.linalg.norm(total))
    return {'paired_output_vector':total.tolist(),'paired_output_norm_rad':norm,
            'candidate_reserve_rad':RESERVE,'beats_candidate_reserve':norm<RESERVE,
            'reserve_margin_rad':RESERVE-norm,'axes':per_axis,
            'pairs_before_norm':True,'source_uniform_verified':False,
            'nonlinear_reset_arithmetic_included':False,'theorem_closed':False}

def run(eigen):
    return {'qualification':'OU3_CARRIED_PAIRED_PHYSICAL_SOURCE_V1',
            'cases':[{'profile':h,**analyze(export(eigen,h))} for h in ('0','wave')],
            'scope':'finite carried p/v/S/a feasibility only; not Q_IMU or source-uniform enclosure'}

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3')); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); out=run(a.eigen); a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps(out,indent=2,sort_keys=True))
