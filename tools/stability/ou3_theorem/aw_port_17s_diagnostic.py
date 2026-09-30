"""17-s carried terminal-AW port diagnostic; non-promoting finite evidence."""
import argparse, hashlib, json, math, subprocess, tempfile
from pathlib import Path
import numpy as np
from .ag_readout_source_diagnostic import REPO, HEADER, instrument

DRIVER=REPO/'tools/stability/aw_port_17s_source.cpp'
TARGET=8.45e-3

def compile_and_run(eigen, heading):
    source=(REPO/HEADER).read_text()
    with tempfile.TemporaryDirectory(prefix='ou3-aw-port-') as d:
        d=Path(d); inc=d/'kalman_ou_iii'; inc.mkdir()
        (inc/HEADER.name).write_text(instrument(source))
        exe=d/'port'
        subprocess.run(['g++','-O2','-std=c++20','-I'+str(d),'-I'+str(REPO/'src'),
                        '-isystem',str(eigen),str(DRIVER),'-o',str(exe)],check=True)
        return json.loads(subprocess.check_output([str(exe),heading],text=True))

def psd_factor(a):
    a=(a+a.T)/2
    w,v=np.linalg.eigh(a); w=np.maximum(w,0)
    return v@np.diag(np.sqrt(w))

def normalized_events(trace):
    out=[]; pending=False; pending_sync=False
    for e in trace['events']:
        k=e['kind']
        if k=='prediction':
            F=np.eye(21); F[:6,:6]=e['F_AG']; F[6:18,6:18]=e['F_LIN']
            F[18:21,18:21]=np.eye(3)*float(e['phi_BA'])
            Q=np.zeros((21,21)); Q[:6,:6]=e['Q_AG']; Q[6:18,6:18]=e['Q_LIN']; Q[18:21,18:21]=e['Q_BA']
            out.append({'kind':'prediction','F':F,'Q':(Q+Q.T)/2,'boundary':False})
            pending=True
        elif k=='sync':
            pending_sync=True
        elif k=='sync_completion':
            before=np.array(e['before'],float); after=np.array(e['after'],float)
            D=(after-before); D=(D+D.T)/2
            out.append({'kind':'prediction','F':np.eye(21),'Q':D,'boundary':pending_sync})
            pending_sync=False; pending=False
        elif k=='correction':
            H=np.array(e['H'],float); R=np.array(e['R'],float); K=np.array(e['K'],float)
            out.append({'kind':'correction','H':H,'R':(R+R.T)/2,'K':K,'boundary':False})
        elif k=='reset':
            d=np.array(e['d'],float).reshape(3)
            x,y,z=d; G=np.eye(21); G[:3,:3]+=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])/2
            out.append({'kind':'reset','G':G,'boundary':False})
    if pending: raise ValueError('unfinished prediction boundary')
    return out

def one_axis(events, axis):
    Y=np.zeros(21); Y[15+axis]=1.0
    action=0.0; samples=[]
    # samples are post-operation backward rows at actual AW sync boundaries.
    for rev,e in enumerate(reversed(events)):
        if e['kind']=='correction':
            L=Y@e['K']; action += float(L@e['R']@L.T)
            Y=Y@(np.eye(21)-e['K']@e['H'])
        elif e['kind']=='prediction':
            action += float(Y@e['Q']@Y.T)
            Y=Y@e['F']
            if e.get('boundary'): samples.append((len(events)-1-rev,Y.copy()))
        else:
            Y=Y@e['G']
    samples.sort()
    # Root nuisance action is diagnostic with literal root covariance, not theorem U_n.
    P0=np.array(trace_global['root_covariance'],float)
    action += float(Y@P0@Y.T)
    hs=[]
    for (idx,Yb),(idx2,_) in zip(samples[:-1],samples[1:]):
        duration=0.0; first=None
        for e in events[idx+1:idx2+1]:
            if e['kind']=='prediction':
                dt=float(e['F'][9,6])
                if dt>0:
                    duration+=dt
                    if first is None: first=e['F']
        if duration<=0 or first is None: continue
        h=np.zeros(3)
        for a in range(3):
            h[a]=(first[6+a,15+a]*Yb[6+a]+first[9+a,15+a]*Yb[9+a]+first[12+a,15+a]*Yb[12+a])/duration
        hs.append(h)
    if not hs: raise ValueError('no complete AW-sync slabs')
    d2=[hs[0]]+[hs[i]-hs[i-1] for i in range(1,len(hs))]+[hs[-1]]
    d21=sum(float(np.linalg.norm(x)) for x in d2)
    root=math.sqrt(max(action,0.0))
    c=d21/root if root>0 else math.inf
    return {'axis':axis,'sync_slabs':len(hs),'reader_action':action,'sqrt_reader_action':root,
            'D2_block_2_1':d21,'C_port_W':c,'four_C_port_W':4*c,
            'target_C_port':TARGET,'ratio_to_target':c/TARGET,'below_target':c<TARGET}

def run(eigen):
    global trace_global
    cases=[]
    for heading in ('0','wave'):
        trace_global=compile_and_run(eigen,heading)
        ev=normalized_events(trace_global)
        axes=[one_axis(ev,a) for a in range(3)]
        cases.append({'profile':heading,'events':len(ev),'axes':axes,
                      'max_C_port_W':max(x['C_port_W'] for x in axes),
                      'all_axes_below_target':all(x['below_target'] for x in axes)})
    return {'qualification':'OU3_AW_PORT_17S_CARRIED_DIAGNOSTIC_V1',
            'window_s':17,'definition':'C_port_W = block-(2,1) norm of D2 h divided by sqrt(causal reader action)',
            'target_C_port':TARGET,'non_promoting':True,'source_uniform_verified':False,
            'theorem_closed':False,'cases':cases}

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3')); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.write_text(json.dumps(run(a.eigen),indent=2,sort_keys=True)+'\n')
