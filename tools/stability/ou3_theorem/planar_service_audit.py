"""Audit complete finite planar operation streams without theorem promotion.

The local S/prediction commutator is evaluated at the true pre-prediction P
with the same literal F,Q,R_S. It is not the complete neighboring scheduler
word defect: nonmagnetic corrections, AW sync, resets and mean feedback remain.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
from .planar_service_stream import records, expand
from .planar_parity import EVEN, ODD
from .planar_service_cell import aw_floor, linked_shift, linked_relative_cell_bound
from .scheduler_s_shift_bound import defect, norm_bound


def symmetric(a):
    return (a+a.T)/2


def reset_matrix(dtheta):
    x,y,z=.5*np.asarray(dtheta)
    G=np.eye(21)
    G[:3,:3]=[[1,-z,y],[z,1,-x],[-y,x,1]]
    return G


def correction_record(a):
    return (expand(a[:225]), a[225:288].reshape(3,21), a[288:297].reshape(3,3),
            a[297:306].reshape(3,3), a[306:369].reshape(21,3),
            a[369:432].reshape(21,3), a[432:435])


def audit(path: Path, stride: int = 1):
    if stride < 1:
        raise ValueError('stride must be positive')
    maxima={};counts=Counter();lastP=None;expected_reset=None;expected_postreset=None
    innovation_min=1e100;cov_min=1e100;sync_active_min=1e100;sync_active_max=-1e100
    extrema={};shift={};samples=[];due_samples=set();pred_samples=set();local_by_sample={}
    def err(name,a):
        maxima[name]=max(maxima.get(name,0.),float(np.max(np.abs(a))))
    for kind,k,a in records(path):
        counts[kind]+=1
        if kind==1:
            P,F,Q,Pminus=(expand(a[i:i+225]) for i in (0,225,450,675))
            if lastP is not None: err('next_prediction_input',P-lastP)
            err('prediction',Pminus-(F@P@F.T+Q));lastP=Pminus
            pred_samples.add(k)
            if (counts[1]-1)%stride==0:
                row={};R=a[900:909].reshape(3,3)
                for name,idx,axes in [('even',EVEN,(0,2)),('odd',ODD,(1,))]:
                    ix=np.ix_(idx,idx);p=symmetric(P[ix]);f=F[ix];q=symmetric(Q[ix])
                    h=np.zeros((len(axes),len(idx)))
                    for j,axis in enumerate(axes):h[j,idx.index(12+axis)]=1
                    r=R[np.ix_(axes,axes)]
                    root=np.linalg.cholesky(p);W=np.linalg.solve(root,np.eye(len(idx)))
                    D,info=linked_shift(p,f,q,h,r)
                    Dw,winfo=linked_shift(p,f,q,h,r,W)
                    err('linked_'+name+'_identity',D-defect(p,f,q,h,r))
                    row[name]={'absolute':info['norm_diagnostic'], 'relative':winfo['norm_diagnostic'],
                        'two_term_bound':norm_bound(p,p,f,q,h,r),
                        'linked_point_bound':linked_relative_cell_bound(p,root,0.,f,q,h,r,W),
                        'linked_one_percent_cell_bound':linked_relative_cell_bound(p,root,.01,f,q,h,r,W),
                        'rank_upper':info['rank_upper']}
                local_by_sample[k]=row
        elif kind in (2,3,4):
            P,H,R,S,K,PCt,res=correction_record(a)
            if lastP is not None:err('correction_input',P-lastP)
            err('innovation',S-(H@P@H.T+R));err('PHt',PCt-P@H.T);err('gain_solve',K@S-PCt)
            innovation_min=min(innovation_min,float(np.linalg.eigvalsh(symmetric(S)).min()))
            expected_reset=P-K@PCt.T-PCt@K.T+K@S@K.T
            if kind==4:due_samples.add(k)
        elif kind==5:
            P=expand(a[:225]);G=reset_matrix(a[225:228])
            if expected_reset is None:raise ValueError('reset without correction')
            err('joseph',P-expected_reset);expected_postreset=G@P@G.T;expected_reset=None
        elif kind==8:
            P=expand(a)
            if expected_postreset is None:raise ValueError('reset output without input')
            err('reset_congruence',P-expected_postreset);lastP=P;expected_postreset=None
        elif kind==6:
            P=expand(a[:225]);target=a[225:234].reshape(3,3);pending=bool(a[234])
            if pending:
                gap=symmetric(target-lastP[15:18,15:18])
                vals=np.linalg.eigvalsh(gap);sync_active_min=min(sync_active_min,float(vals.min()))
                sync_active_max=max(sync_active_max,float(vals.max()));counts['aw_pending']+=1
                counts['aw_fully_active' if vals.min()>0 else 'aw_not_fully_active']+=1
                expected=aw_floor(lastP,target,(15,16,17))
            else:expected=lastP
            err('aw_sync',P-expected);lastP=P
        elif kind==7:
            P=expand(a[:225]);err('sample_output',P-lastP);lastP=P;samples.append(k)
            cov_min=min(cov_min,float(np.linalg.eigvalsh(symmetric(P)).min()))
            for key,v in zip(('tau','elapsed','period','tau_applied','sigma_applied','RS_applied',
                    'tau_target','sigma_target','RS_target','wave_period','acc_var','vertical_acc','AW_pending'),a[280:293]):
                pair=extrema.setdefault(key,[float(v),float(v)]);pair[0]=min(pair[0],float(v));pair[1]=max(pair[1],float(v))
    if expected_reset is not None or expected_postreset is not None:raise ValueError('incomplete final correction')
    if not samples or set(samples)!=pred_samples:raise ValueError('incomplete sample coverage')
    for group,keys in [('all_tested_predictions',set(local_by_sample)),('actually_S_due_tested_predictions',due_samples&set(local_by_sample))]:
        target=shift[group]={}
        for parity in ('even','odd'):
            fields={}
            for field in ('absolute','relative','two_term_bound','linked_point_bound','linked_one_percent_cell_bound'):
                if keys:
                    k=max(keys,key=lambda k:local_by_sample[k][parity][field])
                    fields[field]={'max':local_by_sample[k][parity][field],'sample':k}
            target[parity]=fields
        target['samples']=len(keys)
    return {'qualification':'OU3_PLANAR_SERVICE_OPERATION_AUDIT_V1','result_type':'FINITE DIAGNOSTIC ONLY',
        'stream_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sample_range':[samples[0],samples[-1]],
        'counts':{str(k):v for k,v in counts.items()},'maximum_absolute_operation_residuals':maxima,
        'smallest_covariance_eigenvalue_observed':cov_min,'smallest_innovation_eigenvalue_observed':innovation_min,
        'AW_target_minus_prior_eigenvalue_range':[sync_active_min,sync_active_max],
        'joint_tuner_marginal_ranges_diagnostic_only':extrema,'local_S_prediction_commutator':shift,
        'local_commutator_is_complete_scheduler_word_defect':False,
        'finite_float32_residuals_are_not_uniform_arithmetic_bounds':True,
        'scheduler_phase_cell_forward_invariant':False,'all_time_magnetic_service_verified':False,
        'theorem_closed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stream',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--stride',type=int,default=1)
    a=p.parse_args();out=audit(a.stream,a.stride);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
