"""Actual paired shipping words, not a source-uniform stability certificate.

The metric is fixed below before executing probes; it is never fitted to tapes.
Finite secants are not Jacobian columns. Root variants have different covariance
and frontend ancestry, and must not be assembled into a fictitious 24x24 map.
80/120-digit evaluation checks graph algebra on recorded binary32 operands;
it does not re-execute Eigen or certify its roundoff over a source family.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import open_obligations
from diagnose import bias_contracts

# Predeclared diagnostic units, NOT certified basin radii or physical premises.
SCALES = np.repeat([.1, .01, 1., 1., 10., 1., .1, .1], 3)
METRIC = np.diag(1/SCALES**2)
for j in range(3):
    METRIC[18+j,21+j] = METRIC[21+j,18+j] = 25.
DT = .005
OMEGA = 2*math.pi/8


@dataclass
class Tape:
    tags: np.ndarray
    data: np.ndarray
    family: int
    root_axis: int
    epsilon: float
    amplitude: float
    start: int
    live: int
    direction: int


def load(path: Path) -> Tape:
    meta = path.with_suffix('.meta').read_text().split()
    if len(meta) != 7:
        raise ValueError('incomplete source provenance')
    tags = np.loadtxt(path, delimiter=',', usecols=0, dtype=str)
    data = np.loadtxt(path, delimiter=',', usecols=range(1,38))
    if data.ndim != 2 or not np.isfinite(data[:,:32]).all() or np.isinf(data).any():
        raise ValueError('invalid/nonfinite state tape')
    t = Tape(tags, data, int(meta[0]), int(meta[1]), float(meta[2]),
             float(meta[3]), int(meta[4]), int(meta[5]), int(meta[6]))
    if tags[0] != 'word_start' or tags[-1] != 'word_end':
        raise ValueError('incomplete literal word')
    if data[0,0] != t.start or data[-1,0] != t.start+600:
        raise ValueError('wrong word duration')
    return t


def qmul(a, b):
    """Scalar-first quaternion product, broadcast over leading dimensions."""
    return np.concatenate((a[...,:1]*b[...,:1]-(a[...,1:]*b[...,1:]).sum(axis=-1,keepdims=True),
        a[...,:1]*b[...,1:]+b[...,:1]*a[...,1:]+np.cross(a[...,1:],b[...,1:])),axis=-1)


def effective_quaternion(data):
    # Real-arithmetic extension across pending attitude injection. Actual float
    # injection/normalization/reset changes remain visible in the next event.
    d=data[:,6:9]; theta=np.linalg.norm(d,axis=1); t2=theta**2
    small=theta < .01
    w=np.where(small,1-t2/8+t2*t2/384,np.cos(theta/2))
    k=np.where(small,.5-t2/48+t2*t2/3840,
               np.sin(theta/2)/np.where(theta==0,1,theta))
    corr=np.column_stack((w,k[:,None]*d))
    corr/=np.linalg.norm(corr,axis=1)[:,None]
    q=qmul(corr,data[:,2:6]); q/=np.linalg.norm(q,axis=1)[:,None]
    return q


def errors(tape: Tape):
    """True-minus-estimate joint24 including one-time physical Live primitive."""
    a=tape.data; t=a[:,0]*DT; s=OMEGA*t
    roll=.08*np.sin(s); pitch=.05*np.sin(s+.7)
    qbw=Rotation.from_euler('ZYX',np.column_stack((np.full(len(t),.3),pitch,roll)))
    q=effective_quaternion(a)
    qest=Rotation.from_quat(q[:,[1,2,3,0]])
    eatt=(qbw.inv()*qest.inv()).as_rotvec()
    unit=np.column_stack((.25*np.sin(s),.15*np.cos(s),np.sin(s+.3)))
    p=tape.amplitude*unit
    v=tape.amplitude*OMEGA*np.column_stack((.25*np.cos(s),-.15*np.sin(s),np.cos(s+.3)))
    aw=-tape.amplitude*OMEGA**2*unit
    def primitive(x):
        return tape.amplitude/OMEGA*np.column_stack((-.25*np.cos(x),.15*np.sin(x),-np.cos(x+.3)))
    S=primitive(s)-primitive(np.array([OMEGA*tape.live*DT]))
    b = (.01*np.exp(-t/1200)+.015+.004*np.sin(s/75) if tape.family==0 else
         .04*np.exp(-t/1200)+.008*np.sin(s/75) if tape.family==1 else
         .02+.005*np.sin(s/75))
    beta=b[:,None]*np.array([1,-.5,.25])
    if tape.root_axis>=0:
        beta[:,tape.root_axis]+=tape.epsilon*.1*(np.exp(-t/1200) if tape.family==1 else 1)
    return np.column_stack((eatt,-a[:,9:12],v-a[:,12:15],p-a[:,15:18],
                            S-a[:,18:21],aw-a[:,21:24],beta-a[:,24:27],beta))


def energy(z):
    return np.einsum('...i,ij,...j->...',z,METRIC,z)


def compare(base: Tape, other: Tape):
    if (base.family,base.epsilon,base.amplitude,base.start) != (other.family,other.epsilon,other.amplitude,other.start):
        raise ValueError('mixed source provenance or word start')
    if other.root_axis<0 and base.live!=other.live:
        raise ValueError('same-source perturbation changed the one-time Live origin')
    # Retain mismatches rather than pretending two different hybrid words align.
    same_events=np.array_equal(base.tags,other.tags) and np.array_equal(base.data[:,0],other.data[:,0])
    same_guards=same_events and np.array_equal(base.data[:,1],other.data[:,1])
    zb,zo=errors(base),errors(other)
    d0=zo[0]-zb[0]; dn=zo[-1]-zb[-1]
    v0=float(energy(d0)); vn=float(energy(dn))
    if not v0>0: raise ValueError('zero incremental entry energy')
    report={
        'direction':other.direction,'root_axis':other.root_axis,
        'initial_joint24_increment':d0.tolist(),'final_joint24_increment':dn.tolist(),
        'initial_energy':v0,'final_energy':vn,'unsupplied_secant_ratio':vn/v0,
        'distance_to_one':1-vn/v0,'same_event_word':same_events,'same_mode_word':same_guards,
        'initial_neutral_bias_increment':d0[18:].tolist(),
        'physical_truth_shared':other.root_axis<0,
        'source_root_changed_at_boot':other.root_axis>=0,
        'actual_full_joint24_map_certified':False,
        'supply_subtracted':False,
        'zero_noise_sensor_model_defect_enclosed':False,
        'end_modes':[int(other.data[0,1]),int(other.data[-1,1])],
        'max_attitude_error_rad':float(np.max(np.linalg.norm(zo[:,:3],axis=1))),
    }
    if same_events:
        v=energy(zo-zb); changes=np.diff(v)
        signed=defaultdict(float)
        for tag,delta in zip(other.tags[1:],changes): signed[str(tag)]+=float(delta)
        telescope=math.fsum(signed.values())-(vn-v0)
        if abs(telescope)>1e-10*max(v0,vn,1e-12): raise ValueError('signed ledger fails to telescope')
        report.update(signed_delta_by_operation=dict(signed),
                      telescope_error=telescope,largest_prefix_ratio=float(v.max()/v0),
                      frontend_endpoint_max_abs_increment=float(np.nanmax(np.abs(other.data[:,32:]-base.data[:,32:]))))
    else:
        report['signed_delta_by_operation']=None
        report['prefix_alignment_failure']='different literal event words; endpoint comparison only'
    return report


def mp_energy_of_endpoint(tape: Tape, row: int):
    """High-precision evaluation of the declared chart on actual recorded values."""
    x=list(map(mp.mpf,map(float,tape.data[row])))
    t=x[0]*mp.mpf(1)/200; omega=mp.pi/4; s=omega*t
    def mul(a,b):
        w,u=a[0],a[1:]; v,z=b[0],b[1:]
        cross=[u[1]*z[2]-u[2]*z[1],u[2]*z[0]-u[0]*z[2],u[0]*z[1]-u[1]*z[0]]
        return [w*v-sum(i*j for i,j in zip(u,z))]+[w*z[k]+v*u[k]+cross[k] for k in range(3)]
    def normq(q):
        n=mp.sqrt(sum(v*v for v in q)); return [v/n for v in q]
    def inv(q): return [q[0]]+[-v for v in q[1:]]
    def axis(angle,k):
        a=[mp.cos(angle/2),mp.mpf(0),mp.mpf(0),mp.mpf(0)];a[k+1]=mp.sin(angle/2);return a
    roll=mp.mpf('.08')*mp.sin(s); pitch=mp.mpf('.05')*mp.sin(s+mp.mpf('.7'))
    bw=mul(mul(axis(mp.mpf('.3'),2),axis(pitch,1)),axis(roll,0))
    theta=mp.sqrt(sum(v*v for v in x[6:9])); t2=theta**2
    if theta < mp.mpf('.01'): w=1-t2/8+t2*t2/384;k=mp.mpf('.5')-t2/48+t2*t2/3840
    else: w=mp.cos(theta/2);k=mp.sin(theta/2)/theta
    qe=normq(mul(normq([w]+[k*v for v in x[6:9]]),x[2:6]))
    dq=normq(mul(inv(bw),inv(qe)))
    if dq[0]<0:dq=[-v for v in dq]
    n=mp.sqrt(sum(v*v for v in dq[1:]));angle=2*mp.atan2(n,dq[0])
    att=[v*angle/n for v in dq[1:]] if n else [mp.mpf(0)]*3
    amp=mp.mpf(str(tape.amplitude)); p=[amp*mp.mpf('.25')*mp.sin(s),amp*mp.mpf('.15')*mp.cos(s),amp*mp.sin(s+mp.mpf('.3'))]
    vel=[amp*omega*mp.mpf('.25')*mp.cos(s),-amp*omega*mp.mpf('.15')*mp.sin(s),amp*omega*mp.cos(s+mp.mpf('.3'))]
    def prim(u):return [-amp/omega*mp.mpf('.25')*mp.cos(u),amp/omega*mp.mpf('.15')*mp.sin(u),-amp/omega*mp.cos(u+mp.mpf('.3'))]
    S=[i-j for i,j in zip(prim(s),prim(omega*tape.live/200))]
    aw=[-omega**2*v for v in p]
    b=(mp.mpf('.01')*mp.exp(-t/1200)+mp.mpf('.015')+mp.mpf('.004')*mp.sin(s/75) if tape.family==0 else
       mp.mpf('.04')*mp.exp(-t/1200)+mp.mpf('.008')*mp.sin(s/75) if tape.family==1 else
       mp.mpf('.02')+mp.mpf('.005')*mp.sin(s/75))
    beta=[b,-b/2,b/4]
    if tape.root_axis>=0:beta[tape.root_axis]+=mp.mpf(str(tape.epsilon))/10*(mp.exp(-t/1200) if tape.family==1 else 1)
    return mp.matrix(att+[-v for v in x[9:12]]+
        [v-x[12+j] for j,v in enumerate(vel)]+[v-x[15+j] for j,v in enumerate(p)]+
        [v-x[18+j] for j,v in enumerate(S)]+[v-x[21+j] for j,v in enumerate(aw)]+
        [v-x[24+j] for j,v in enumerate(beta)]+beta)


def high_precision_pair(base,other):
    out={}
    for dps in (80,120):
        with mp.workdps(dps):
            m=mp.matrix(METRIC.tolist())
            d0=mp_energy_of_endpoint(other,0)-mp_energy_of_endpoint(base,0)
            dn=mp_energy_of_endpoint(other,-1)-mp_energy_of_endpoint(base,-1)
            out[str(dps)]=mp.nstr((dn.T*m*dn)[0]/(d0.T*m*d0)[0],dps-8)
    out['nonlinear_shipping_execution_precision']='binary32, not 80/120 digits'
    return out


def solve_graph(base_path: Path, other_path: Path, dps=80):
    """No inverse reconstruction from covariance: read actual N,S,K,r."""
    def rows(path):
        with path.open() as f:
            for line in f:
                cells=line.strip().split(','); yield cells[:3],list(map(float,cells[3:]))
    def unpack(v):
        def mat(a,n,m):return mp.matrix([v[a+i*m:a+(i+1)*m] for i in range(n)])
        if len(v)!=621:raise ValueError(f'wrong actual solve record: {len(v)}')
        return (mat(0,3,1),mat(3,21,3),mat(66,3,3),mat(75,21,3),
                mat(138,21,21),mat(579,21,1),mat(600,21,1))
    maxima=defaultdict(float); by_tag=defaultdict(lambda:defaultdict(float)); counts=Counter()
    with mp.workdps(dps):
        from itertools import zip_longest
        for a,b in zip_longest(rows(base_path),rows(other_path)):
            if a is None or b is None or a[0]!=b[0]:
                return {'aligned':False,'reason':'solve occurrence/guard mismatch',
                        'aligned_prefix_pairs':sum(counts.values()),
                        'first_mismatch_base':a[0] if a else None,
                        'first_mismatch_other':b[0] if b else None, 'source_uniform':False}
            tag=a[0][0];counts[tag]+=1
            r0,N0,S0,K0,P0,x0,y0=unpack(a[1]);r1,N1,S1,K1,P1,x1,y1=unpack(b[1])
            # Eigen LDLT consumes a triangle. The full-S ideal graph discrepancy
            # and asymmetric/numerical differences are recorded, not erased.
            q0=mp.lu_solve(S0,r0);q1=mp.lu_solve(S1,r1);dq=q1-q0
            termS=(S1-S0)*q0;termN=(N1-N0)*q0
            ideal=N1*dq+termN; recorded_Kr=K1*r1-K0*r0
            actual=(y1-x1)-(y0-x0)
            quantities={
                'nominal_residual_norm':mp.norm(r0),
                'delta_S_times_q0_norm':mp.norm(termS),
                'delta_N_times_q0_norm':mp.norm(termN),
                'descriptor_residual_norm':mp.norm(S1*dq+termS-(r1-r0)),
                'correction_identity_residual_norm':mp.norm(ideal-(N1*q1-N0*q0)),
                'recorded_Kr_minus_ideal_Nq_increment_norm':mp.norm(recorded_Kr-ideal),
                'actual_mean_minus_recorded_Kr_increment_norm':mp.norm(actual-recorded_Kr),
                'actual_mean_minus_ideal_Nq_increment_norm':mp.norm(actual-ideal),
                'frozen_gain_omitted_correction_norm':mp.norm(ideal-N0*mp.lu_solve(S0,r1-r0)),
                'delta_covariance_frobenius_norm':mp.norm(P1-P0),
                'innovation_asymmetry_norm':mp.norm(S1-S1.T),
            }
            for key,value in quantities.items():
                maxima[key]=max(maxima[key],float(value));by_tag[tag][key]=max(by_tag[tag][key],float(value))
    return {'aligned':True,'decimal_digits':dps,'counts':dict(counts),'maxima':dict(maxima),
            'by_operation':{k:dict(v) for k,v in by_tag.items()},
            'actual_K_and_N_from_shipping':True,'source_uniform':False,
            'deployment_roundoff_enclosed':False,
            'norm_units':'raw mixed state/residual coordinates; not a physical performance bound',
            'qualification':'exact solve identities evaluated on actual recorded operands; measured defects only'}


def analytic_probe_envelope(amplitude=1.2,epsilon=.001):
    # Triangle/component bounds hold for ALL t for these particular analytic
    # probes, not just sampled points; they do not cover all COMPLETE-BRMM.
    n=math.sqrt(1+.25**2+.15**2);p=amplitude*n;v=OMEGA*p;a=OMEGA*v
    return {'method':'analytic all-time component/triangle bounds for this probe subfamily',
            'p_norm_upper_m':p,'v_norm_upper_mps':v,'a_norm_upper_mps2':a,
            'primitive_difference_upper_ms':2*p/OMEGA,
            'omega_body_norm_upper_radps':OMEGA*math.hypot(.08,.05),
            'harmonic_frequency_hz':1/8,
            'kinematic_caps_satisfied':p<=8.1 and v<=5.5 and a<=8.8 and 2*p/OMEGA<=1100,
            'BIAS0_component_bound':.01+.015+.004+abs(epsilon)*.1,
            'BIAS1_component_bound':.04+.008+abs(epsilon)*.1,
            'BIAS2_component_bound':.02+.005+abs(epsilon)*.1,
            'BIAS2_derivative_bound':.005*OMEGA/75,
            'BIAS0_driver':'GM root .01, tau1200; turn-on .015+root increment; strain .004 sin(2*pi*t/600)',
            'BIAS1_driver':'root .04+root increment, tau1200; sinusoid amplitude .008, period600',
            'BIAS2_driver':'phi_true=1, w=beta_next-beta; constant boot root plus slow sinusoid',
            'full_BRMM_PE_sensor_startup_admission':False,'full_source_cover':False}


def run(directory:Path, with_solves=False):
    result={'qualification':'ACTUAL_PAIRED_SHIPPING_FINITE_WORD_DIAGNOSTIC_ONLY',
        'status':open_obligations(),'source_uniform':False,
        'metric':{'scales':SCALES.tolist(),'normalized_bias_cross':.25,
                  'chosen_before_data':True,'fitted_to_replay':False,
                  'coercivity_eigenvalues':[float(v) for v in np.linalg.eigvalsh(METRIC)[[0,-1]]]},
        'supply_budget_certified':False,'common_storage_master_solved':False,
        'source':analytic_probe_envelope(load(directory/'H18_-1.csv').amplitude,
                                         load(directory/'H18_-1.csv').epsilon),
        'bias_contracts':bias_contracts(),
        'modes':{},'trace_sha256':{}}
    for mode in ('H18','A21','edge'):
        basepath=directory/f'{mode}_-1.csv';base=load(basepath)
        expected=(0,0) if mode=='H18' else (1,1) if mode=='A21' else (0,1)
        if tuple(base.data[[0,-1],1].astype(int))!=expected:
            raise ValueError(f'{mode}: missing natural mode/transition')
        records=[]
        for suffix in [str(j) for j in range(21)]+[f'beta{j}' for j in range(3)]:
            path=directory/f'{mode}_{suffix}.csv';other=load(path)
            rec=compare(base,other); rec['file']=path.name; records.append(rec)
        worst=max(records,key=lambda r:r['unsupplied_secant_ratio'])
        motion=max((r for r in records if 0<=r['direction']<18),key=lambda r:r['unsupplied_secant_ratio'])
        otherpath=directory/worst['file'];other=load(otherpath)
        worst['high_precision_chart_evaluation']=high_precision_pair(base,other)
        info={'baseline_events':dict(Counter(base.tags)),'worst_tested_direction':worst,
              'worst_tested_motion_direction':motion,'directions':records,
              'basis_secants_not_a_Jacobian':True,'hybrid_branch_mismatches':sum(not r['same_mode_word'] for r in records),
              'every_prefix_chart_retention_certified':False}
        if with_solves:
            # Detailed direction is selected by observation, not used to fit M.
            pairs=[]
            for path in directory.glob(f'{mode}_*.solves.csv'):
                if path.name==f'{mode}_-1.solves.csv':continue
                pairs.append((path,solve_graph(directory/f'{mode}_-1.solves.csv',path)))
            info['actual_solve_graphs']={p.name:r for p,r in pairs}
            selected=directory/(worst['file'][:-4]+'.solves.csv')
            if selected.exists():
                hi=solve_graph(directory/f'{mode}_-1.solves.csv',selected,dps=120)
                info['selected_graph_120_digits']=hi
                low=info['actual_solve_graphs'].get(selected.name)
                if low and low['aligned'] and hi['aligned']:
                    info['graph_80_120_max_defect_difference']=abs(
                        hi['maxima']['actual_mean_minus_ideal_Nq_increment_norm']-
                        low['maxima']['actual_mean_minus_ideal_Nq_increment_norm'])
        result['modes'][mode]=info
        for path in (basepath,otherpath):result['trace_sha256'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    return result



def control_equivalence(observed: Path, control: Path):
    """Strict same-compiler comparison against an unmodified shipping header."""
    pa=json.loads((observed.parent/'provenance.json').read_text())
    pb=json.loads((control.parent/'provenance.json').read_text())
    if pa.get('instrumented') is not True or pb.get('instrumented') is not False:
        raise ValueError('control must use the uninstrumented build')
    if pa['shipping_header_sha256'] != pb['shipping_header_sha256']:
        raise ValueError('control and observed shipping sources differ')
    if pa['compile_command'][:4]+pa['compile_command'][5:9] != pb['compile_command'][:4]+pb['compile_command'][5:9]:
        raise ValueError('control compiler/options differ')
    paths=sorted(p for p in control.glob('*.csv') if not p.name.endswith('.solves.csv'))
    if len(paths)!=75: raise ValueError('control must contain all 75 words')
    count=0
    for path in paths:
        a=load(observed/path.name);b=load(path)
        keep=np.isin(a.tags,['word_start','word_end','wrapper_imu','wrapper_mag'])
        if not np.array_equal(a.tags[keep],b.tags) or not np.array_equal(a.data[keep].view(np.uint64),b.data.view(np.uint64)):
            raise ValueError(f'control differs from observed shipping execution: {path.name}')
        count+=len(b.tags)
    return {'bitwise_equal_recorded_wrapper_states':True,'words':len(paths),'rows':count,
            'scope':'same compiler/options; not deployment arithmetic qualification'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--solves',action='store_true')
    ap.add_argument('--control-directory',type=Path)
    a=ap.parse_args();report=run(a.directory,a.solves)
    if a.control_directory:
        report['uninstrumented_control']=control_equivalence(a.directory,a.control_directory)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:{'ratio':v['worst_tested_direction']['unsupplied_secant_ratio'],
                       'direction':v['worst_tested_direction']['file'],
                       'motion_ratio':v['worst_tested_motion_direction']['unsupplied_secant_ratio'],
                       'branch_mismatches':v['hybrid_branch_mismatches']} for k,v in report['modes'].items()},indent=2))

if __name__=='__main__':main()
