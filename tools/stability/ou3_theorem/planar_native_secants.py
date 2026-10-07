"""Finite central secants of the FULL shipping execution, not derivative bounds.

The native probe forks the entire inherited execution at sample 40000, perturbs
all 21 mean and 231 symmetric covariance directions, then continues the original
shipping wrapper on the same delivered record. This includes nonlinear gains,
Joseph/reset, projections, frontend, adaptation and clocks. Perturbed root
states are diagnostic neighbourhood points; reachability is NOT asserted.
"""
from __future__ import annotations
import argparse,hashlib,json,struct,subprocess,tempfile
from pathlib import Path
import numpy as np
from .moving_compatibility_diagnostic import REPO
from .compatibility_quotient import whitened_line
from .kernel_restricted_action import householder_complement
from .planar_compatibility_quotient_mean import compatibility_line
PROBE=Path(__file__).with_name('planar_native_secant_probe.cpp')


def read(path):
    data=Path(path).read_bytes()
    if len(data)!=2*466*4+12: raise ValueError('incomplete native secant record')
    a=np.frombuffer(data[:2*466*4],dtype='<f4').astype(float).reshape(2,466)
    if not np.isfinite(a).all(): raise ValueError('nonfinite secant state')
    return a,struct.unpack('<QI',data[-12:])


def difference(state,reference):
    x=state[:21]-reference[:21]
    q=state[21:25]; q=q/np.linalg.norm(q)
    r=reference[21:25]; r=r/np.linalg.norm(r)
    w=q@r; v=r[0]*q[1:]-q[0]*r[1:]-np.cross(q[1:],r[1:])
    if w<0: w=-w; v=-v
    n=np.linalg.norm(v)
    x[:3]=(2*np.arctan2(n,w)/n)*v if n>0 else 0.
    return x


def svec(A):
    return np.array([A[i,j]*(1. if i==j else np.sqrt(2.)) for i in range(21) for j in range(i,21)])


def analyze(directory,epsilon):
    if not np.isfinite(epsilon) or not 0<epsilon<.1: raise ValueError('finite diagnostic epsilon in (0,.1) required')
    directory=Path(directory); nominal,forcing=read(directory/'nominal.bin')
    control,control_forcing=read(directory/'-1_-1.bin')
    if not np.array_equal(nominal,control) or forcing!=control_forcing:
        raise ValueError('unperturbed fork differs from parent execution')
    P0=nominal[0,25:].reshape(21,21); PN=nominal[1,25:].reshape(21,21)
    L0,u0,_,_,_=whitened_line(P0,compatibility_line(40000))
    LN,uN,_,_,_=whitened_line(PN,compatibility_line(44000))
    W0=np.linalg.solve(L0,np.eye(21)); WN=np.linalg.solve(LN,np.eye(21))
    U0=householder_complement(u0); UN=householder_complement(uN)
    D=np.zeros((252,252)); Dc=np.zeros_like(D); root_defect=0.; hash_all=hashlib.sha256()
    def coords(s,ref,W):
        return np.r_[W@difference(s,ref),svec(W@(s[25:].reshape(21,21)-ref[25:].reshape(21,21))@W.T)]
    for j in range(252):
        plus,pf=read(directory/f'{j}_1.bin'); minus,mf=read(directory/f'{j}_-1.bin')
        if pf!=forcing or mf!=forcing: raise ValueError('generated forcing/clock or accepted cadence changed')
        for sign in (-1,1): hash_all.update((directory/f'{j}_{sign}.bin').read_bytes())
        D[:,j]=(coords(plus[1],nominal[1],WN)-coords(minus[1],nominal[1],WN))/(2*epsilon)
        Dc[:,j]=(coords(plus[1],nominal[1],W0)-coords(minus[1],nominal[1],W0))/(2*epsilon)
        initial=(coords(plus[0],nominal[0],W0)-coords(minus[0],nominal[0],W0))/(2*epsilon)
        initial[j]-=1; root_defect=max(root_defect,float(np.linalg.norm(initial)))
    A=D[:21,:21]; B=D[21:,:21]; C=D[:21,21:]; R=D[21:,21:]
    norm=lambda x:float(np.linalg.norm(x,2))
    rhoP=norm(R); rhoQ=norm(UN.T@A@U0); b=norm(B@U0); c=norm(UN.T@C)
    # A tube about a repeating root anchor needs a COMMON norm plus the
    # actual root-to-end drift; endpoint-varying metrics alone do not suffice.
    commonP=norm(Dc[21:,21:]); commonQ=norm(U0.T@Dc[:21,:21]@U0)
    commonb=norm(Dc[21:,:21]@U0); commonc=norm(U0.T@Dc[:21,21:])
    drift=coords(nominal[1],nominal[0],W0)
    result={'epsilon':epsilon,'mean_dimension':21,'covariance_dimension':231,
        'full_mean_secant_gain':norm(A),'quotient_mean_secant_gain':rhoQ,
        'covariance_secant_gain':rhoP,'b_quotient_to_covariance_secant':b,'c_covariance_to_quotient_secant':c,
        'gauge_to_transverse_secant':float(np.linalg.norm(UN.T@A@u0)),
        'gauge_to_covariance_secant':float(np.linalg.norm(B@u0)),
        'comparison_determinant_diagnostic':(1-rhoP)*(1-rhoQ)-b*c,
        'realized_initial_secant_identity_defect_max':root_defect,
        'unperturbed_fork_bitwise_equal':True,'forcing_clock_hash_equal_for_all_perturbations':True,
        'accepted_magnetic_corrections_each_20s_word':forcing[1],
        'native_secant_bytes_sha256':hash_all.hexdigest()}
    result['common_root_metric']={'rho_P_secant':commonP,'rho_Q_secant':commonQ,
        'b_Q_secant':commonb,'c_Q_secant':commonc,
        'determinant_diagnostic':(1-commonP)*(1-commonQ)-commonb*commonc,
        'q_P_point_root_to_end_drift':float(np.linalg.norm(drift[21:])),
        'q_Q_point_root_to_end_drift':float(np.linalg.norm(U0.T@drift[:21])),
        'gauge_drift_point':float(u0@drift[:21]),
        'C_Q_secant':float(np.linalg.norm(U0.T@Dc[:21,:21]@u0)),
        'C_P_secant':float(np.linalg.norm(Dc[21:,:21]@u0)),
        'future_forcing_drift_bounded':False,'uniform_affine_charges_certified':False}
    determinant=(1-commonP)*(1-commonQ)-commonb*commonc
    if commonP<1 and commonQ<1 and determinant>0:
        rp=((1-commonQ)*np.linalg.norm(drift[21:])+commonb*np.linalg.norm(U0.T@drift[:21]))/determinant
        result['common_root_metric']['optimistic_zero_gauge_P_radius_diagnostic']=float(rp)
        result['common_root_metric']['relative_SPD_ball_radius_must_be_below']=1.
        result['common_root_metric']['fixed_root_ball_feasibility_failure_class']='D_SUFFICIENT_BOUND_FAILURE' if rp>=1 else None
        result['common_root_metric']['zero_gauge_role']='optimistic feasibility lower charge only; physical gauge forcing is not discarded from any theorem'
    return result,D


def run(eigen,output,cxx='g++',epsilons=(.01,.005)):
    if not (Path(eigen)/'Eigen/Dense').is_file(): raise ValueError('Eigen/Dense missing')
    with tempfile.TemporaryDirectory(prefix='ou3-native-secants-') as tmp:
        tmp=Path(tmp); exe=tmp/'probe'
        cmd=[cxx,'-std=c++20','-O1','-DEIGEN_UNROLLING_LIMIT=0','-DEIGEN_NON_ARDUINO','-I'+str(REPO/'src'),'-I'+str(eigen),str(PROBE),'-o',str(exe)]
        subprocess.run(cmd,check=True,timeout=240)
        results=[]; matrices=[]
        for i,eps in enumerate(epsilons):
            directory=tmp/str(i); directory.mkdir()
            subprocess.run([str(exe),str(eps),str(directory)],check=True,timeout=1800)
            result,D=analyze(directory,eps); results.append(result); matrices.append(D)
        out=report(results,matrices)
    Path(output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    return out


def shipping_digest():
    source=hashlib.sha256()
    for p in sorted((REPO/'src').rglob('*')):
        if p.is_file():source.update(str(p.relative_to(REPO)).encode()+b'\0'+hashlib.sha256(p.read_bytes()).digest())
    return source.hexdigest()


def verify_report(out):
    if out.get('qualification')!='OU3_PLANAR_FULL_NATIVE_SECANTS_V1' or out.get('result_type')!='FINITE DIAGNOSTIC ONLY':
        raise ValueError('native secants must remain finite diagnostics')
    for key in ('shipping_behavior_changed','complete_nonlinear_derivative_certified','uniform_cross_gains_certified',
                'joint_cell_forward_invariant','all_time_magnetic_service_verified','theorem_closed'):
        if out.get(key) is not False: raise ValueError('missing qualification or promotion: '+key)
    for key,path in (('probe_sha256',PROBE),('analysis_driver_sha256',Path(__file__)),
                     ('shipping_header_sha256',REPO/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h')):
        if out.get(key)!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('source binding changed: '+key)
    if out.get('shipping_source_tree_sha256')!=shipping_digest():raise ValueError('shipping source tree changed')
    if out.get('word_samples')!=[40000,44000] or len(out.get('runs',[]))!=2:raise ValueError('missing complete secant/refinement words')
    for row,epsilon in zip(out['runs'],(.01,.005)):
        if row.get('epsilon')!=epsilon or row.get('mean_dimension')!=21 or row.get('covariance_dimension')!=231:
            raise ValueError('wrong secant coordinates')
        if row.get('unperturbed_fork_bitwise_equal') is not True or row.get('forcing_clock_hash_equal_for_all_perturbations') is not True:
            raise ValueError('native controls failed')
    return True


def report(results,matrices):
    return {'qualification':'OU3_PLANAR_FULL_NATIVE_SECANTS_V1','result_type':'FINITE DIAGNOSTIC ONLY',
        'word_samples':[40000,44000],'runs':results,
        'profile':{'sigma_a':.2,'sigma_g':.00135,'sigma_m':.8,'tau_scaled_S_cadence':True,
                   'scope':'inherited planar probe configuration of shipping wrapper; not the distinct AtomS3R sketch profile'},
        'step_halving_matrix_difference_2norm':float(np.linalg.norm(matrices[-1]-matrices[0],2)) if len(matrices)>1 else None,
        'operator_scope':'central secants of complete binary32 wrapper, not a classical derivative or uniform bound',
        'structures_preserved':'complete inherited wrapper state, physical sensor record, all 21 states and full covariance, gains, Joseph, resets, projections, frontend/tuner and clocks',
        'relaxations_introduced':'off-orbit root perturbations are diagnostic neighbourhood points, not certified reachable histories; finite differences are not derivative enclosures',
        'probe_sha256':hashlib.sha256(PROBE.read_bytes()).hexdigest(),
        'analysis_driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'shipping_source_tree_sha256':shipping_digest(),
        'shipping_header_sha256':hashlib.sha256((REPO/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_bytes()).hexdigest(),
        'shipping_behavior_changed':False,'complete_nonlinear_derivative_certified':False,
        'uniform_cross_gains_certified':False,'joint_cell_forward_invariant':False,
        'all_time_magnetic_service_verified':False,'theorem_closed':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3'));p.add_argument('--cxx',default='g++');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); print(json.dumps(run(a.eigen,a.output,a.cxx),sort_keys=True))
if __name__=='__main__':main()
