"""Finite frozen-coefficient covariance sensitivity, with every AW floor retained.

This is NOT the derivative of the complete closed-loop estimator: mean,
frontend and generated-coefficient variations are explicit open input ports.
It tests covariance-cell feasibility before any interval theorem is attempted.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import correction_record,reset_matrix,symmetric
from .planar_parity import EVEN,ODD


def symmetric_basis(n):
    out=[]
    for i in range(n):
        x=np.zeros((n,n));x[i,i]=1;out.append(x)
        for j in range(i+1,n):
            x=np.zeros((n,n));x[i,j]=x[j,i]=2**-.5;out.append(x)
    return np.asarray(out)


def compressed_maps(path,start,end):
    maps=[];T=np.eye(21);root=None;terminal=None;pre_sync=None;active_min=1e100;counts={}
    for kind,k,a in records(path):
        if k<=start:continue
        if k>end:break
        counts[kind]=counts.get(kind,0)+1
        if kind==1:
            if root is None:root=expand(a[:225])
            T=expand(a[225:450])@T;pre_sync=expand(a[675:900])
        elif kind in (2,3,4):
            P,H,R,S,K,PCt,res=correction_record(a)
            T=(np.eye(21)-K@H)@T
        elif kind==5:T=reset_matrix(a[225:228])@T
        elif kind==6 and bool(a[234]):
            target=a[225:234].reshape(3,3)
            floor_gap=np.linalg.eigvalsh(symmetric(target-pre_sync[15:18,15:18])).min()
            if floor_gap<=0:raise ValueError('active-face derivative not applicable')
            active_min=min(active_min,float(floor_gap));maps.append((T,True));T=np.eye(21)
        elif kind==7:terminal=expand(a[:225])
    if root is None or terminal is None:raise ValueError('empty word')
    maps.append((T,False));return root,terminal,maps,counts,active_min


def apply_maps(batch,maps,idx):
    batch=np.array(batch,copy=True)
    aw=[j for j,i in enumerate(idx) if 15<=i<18]
    for A,floor in maps:
        A=A[np.ix_(idx,idx)]
        batch=A@batch@A.T
        if floor:
            for i in aw:
                for j in aw:batch[:,i,j]=0
    return batch


def diagnostic(path,start,end):
    P0,PN,maps,counts,active_min=compressed_maps(path,start,end)
    out={}
    for name,idx in [('even',EVEN),('odd',ODD)]:
        p0=symmetric(P0[np.ix_(idx,idx)]);pn=symmetric(PN[np.ix_(idx,idx)])
        L0=np.linalg.cholesky(p0);LN=np.linalg.cholesky(pn);WN=np.linalg.inv(LN)
        basis=symmetric_basis(len(idx));d0=L0@basis@L0.T
        dn=apply_maps(d0,maps,idx);dn=WN@dn@WN.T
        # Coordinates in orthonormal symmetric-matrix bases preserve Frobenius.
        op=np.einsum('aij,bij->ab',basis,dn,optimize=True)
        singular=np.linalg.svd(op,compute_uv=False)
        out[name]={'symmetric_dimension':len(basis),'relative_Frobenius_gain':float(singular[0]),
            'leading_singular_values':singular[:8].tolist(),
            'same_root_end_covariance_used':False}
    return {'qualification':'OU3_PLANAR_COVARIANCE_FRECHET_DIAGNOSTIC_V1',
        'result_type':'FINITE DIAGNOSTIC ONLY','root_sample':start,'end_sample':end,
        'stream_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        'retained_operations':{str(k):v for k,v in counts.items()},'compressed_congruence_segments':len(maps),
        'observed_active_AW_floor_gap_min':active_min,'parity_blocks':out,
        'all_AW_sync_derivatives_retained':True,'generated_mean_coefficient_variations_included':False,
        'rounding_certified':False,'scheduler_branch_coverage_certified':False,
        'forward_invariant_cell_certified':False,'all_time_magnetic_service_verified':False,
        'theorem_closed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stream',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--start',type=int,default=40000)
    p.add_argument('--end',type=int,default=44000);a=p.parse_args()
    a.output.write_text(json.dumps(diagnostic(a.stream,a.start,a.end),indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
