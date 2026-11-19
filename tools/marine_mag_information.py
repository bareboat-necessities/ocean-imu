#!/usr/bin/env python3
"""Local profiled magnetic information, independent of the navigation estimator."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation
from sim_dataset import ROOT
from marine_mag_study import OUT, write

BASIS=np.array([[0,1,0],[1,0,0],[0,0,-1.]])
FIELD=52*np.array([np.cos(np.deg2rad(66.5))*np.sin(np.deg2rad(-12.6)),
                   np.cos(np.deg2rad(66.5))*np.cos(np.deg2rad(-12.6)),-np.sin(np.deg2rad(66.5))])
# Scale-free coordinates, Frobenius-orthonormal basis.
I=np.eye(3)
DIAG=np.array([np.diag(v) for v in null_space(np.ones((1,3))).T])
SYM=np.array([(np.outer(I[i],I[j])+np.outer(I[j],I[i]))/np.sqrt(2) for i,j in [(0,1),(0,2),(1,2)]])
SKEW=np.array([(np.outer(I[i],I[j])-np.outer(I[j],I[i]))/np.sqrt(2) for i,j in [(0,1),(0,2),(1,2)]])
NAMES=['bx','by','bz','diag1','diag2','symXY','symXZ','symYZ','rotXY','rotXZ','rotYZ','iso']
MODES=[('hard',3),('diagonal',5),('symmetric',8),('affine',11),('affine_free_scale',12)]


def cross(v):
    x,y,z=v.T
    zero=np.zeros_like(x)
    return np.stack([zero,-z,y,z,zero,-x,-y,x,zero],axis=1).reshape(-1,3,3)


def design(c):
    f=(c@FIELD)@BASIS.T
    u=BASIS@(c@c[0].T)@BASIS.T
    # Integral of the Jacobian of gyro transport for one common body-fixed bias.
    bt=np.swapaxes(u,1,2)
    integral=np.zeros_like(bt)
    integral[1:]=np.cumsum((bt[1:]+bt[:-1])*(0.04/2),axis=0)
    bg=-0.001*cross(f)@u@integral
    e=np.concatenate((DIAG,SYM,SKEW,[I/np.sqrt(3)]))
    gains=np.einsum('aij,tj->tia',e,f)*0.01
    return np.concatenate((u,bg,np.broadcast_to(I,u.shape),gains),axis=2)


def profile(j):
    j=j.reshape(-1,j.shape[-1])/0.8
    n=j[:,:6]
    u,s,_=np.linalg.svd(n,full_matrices=False)
    q=u[:,s>s[0]*1e-11]
    return j[:,6:]-q@(q.T@j[:,6:])


def score(z):
    rows=[]
    for name,p in MODES:
        _,s,v=np.linalg.svd(z[:,:p],full_matrices=False)
        tol=max(1e-8,s[0]*1e-10)
        rank=int(np.sum(s>tol))
        weak=v[-1]
        rows.append(dict(model=name,rank=rank,parameters=p,practical=int(np.sum(s>=2)),
                         smax=s[0],smin=s[-1],condition=s[0]/s[-1] if rank==p else 'inf',
                         weak=' '.join(f'{k}={w:.3f}' for k,w in zip(NAMES,weak) if abs(w)>.1),
                         singular=' '.join(f'{x:.7g}' for x in s)))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=OUT)
    args=p.parse_args()
    rows=[]
    cache=ROOT/'.sim-data/cases'
    for case in ['following','head']:
        for height in ['0.270','1.500','8.500']:
            path=next((cache/case).glob(f'wave_data_jonswap_H{height}_*.csv'))
            df=pd.read_csv(path,usecols=['time','yaw_deg','q_wb_zu_x','q_wb_zu_y','q_wb_zu_z','q_wb_zu_w'])
            data=df.iloc[7::8] # first magnetic tick at the eighth 200 Hz sample
            t=data.time.to_numpy()
            actual=Rotation.from_quat(data[['q_wb_zu_x','q_wb_zu_y','q_wb_zu_z','q_wb_zu_w']]).as_matrix()
            heading=data.yaw_deg.to_numpy()
            mean=0 if case=='following' else 180
            for steering in [0,1,3,5,10,'turn']:
                amp=3 if steering=='turn' else steering
                angle=mean+amp*np.sin(2*np.pi*t/20)
                if steering=='turn':
                    angle+=15*(np.tanh((t-550)/12)-np.tanh((t-630)/12))
                for waves in ([True,False] if height=='0.270' else [True]):
                    if waves:
                        c=actual if steering==3 else actual@Rotation.from_euler('z',(heading-angle)[:,None],degrees=True).as_matrix()
                    else:
                        c=Rotation.from_euler('z',-angle[:,None],degrees=True).as_matrix()
                    for window in [30,120,300,1200]:
                        # Use block nuisance resets to limit gyro drift and avoid
                        # counting correlated gyro noise as 25 independent fixes/s.
                        chunks=[]
                        end=window*25
                        for start in range(0,end,500):
                            j=design(c[start:min(start+500,end)])
                            chunks.append(profile(j[::10]))
                        z=np.concatenate(chunks)
                        for row in score(z):
                            rows.append(dict(case=case,height=height,steering=steering,waves=int(waves),
                                             window=window,**row))
    write(args.out/'information.csv',rows)


if __name__=='__main__':
    main()
