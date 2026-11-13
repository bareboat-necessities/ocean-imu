#!/usr/bin/env python3
"""Offline relinearized reference; never feeds a result into the live estimator.

Truth columns are read ONLY in score(), after fitting. Latent CARMA motion,
per-node attitudes and a constant gyro bias are nuisance variables. Samples
are disjoint windows for this diagnostic, not extra live measurements.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg,sparse
from scipy.optimize._numdiff import approx_derivative,group_columns
from types import SimpleNamespace
from scipy.sparse.linalg import spsolve
from scipy.spatial.transform import Rotation as Rot

def whiten(cov):
    cov=(cov+cov.T)/2
    d,U=np.linalg.eigh(cov)
    floor=max(d.max()*1e-13,1e-18)
    return (U/np.sqrt(np.maximum(d,floor))).T,dict(min_eigenvalue=float(d.min()),floor=float(floor),floored=int((d<floor).sum()))

BIAS_PRIOR_SD=.008*9.80665/np.sqrt(3)

def run(prefix,prior_scale=1.,max_nfev=30,bias_prior_sd=BIAS_PRIOR_SD):
    data=pd.read_csv(str(prefix)+'_window.csv');model=json.loads(Path(str(prefix)+'_model.json').read_text())
    # Use actual magnetic timestamps, ten packets per node (~0.4 s).
    ix=np.flatnonzero(data.mag_tick.to_numpy())[::10];d=data.iloc[ix].reset_index(drop=True);n=len(d)
    assert n>40
    dt=np.diff(d.t.to_numpy());assert np.max(abs(dt-dt[0]))<1e-6
    h=dt[0];A=np.kron(np.eye(3),model['A6']);W=np.array(model['W18']);W=(W+W.T)/2
    P=linalg.solve_continuous_lyapunov(A,-W);scale=np.sqrt(np.maximum(P.diagonal(),1e-15))
    As=(A*scale[None,:])/scale[:,None];Ws=(W/scale[None,:])/scale[:,None]
    F=linalg.expm(As*h)
    E=linalg.expm(np.block([[As,Ws],[np.zeros_like(As),-As.T]])*h)
    Q=E[:18,18:]@F.T
    Lq,qdiag=whiten(Q*prior_scale);Lp,pdiag=whiten((P/scale[None,:]/scale[:,None])*prior_scale)
    field=np.array(model['field'])[0];g=np.array([0,0,9.80665])
    acc=d[['ax','ay','az']].to_numpy();mag=d[['mx','my','mz']].to_numpy()
    q0=Rot.from_quat(d[['qx','qy','qz','qw']].to_numpy());m=int(ix[1]-ix[0]);assert np.all(np.diff(ix)==m)
    raw=data[['gx','gy','gz']].to_numpy();rates=np.stack([.5*(raw[ix[k]:ix[k+1]]+raw[ix[k]+1:ix[k+1]+1]) for k in range(n-1)])
    hgyro=.005;sdgyro=model['gyro_sigma']*np.sqrt(m)*hgyro
    # The conservative shared-endpoint gyro-noise bound uses this full value,
    # not sqrt(1/2) times it. No sample noise is counted as independent twice.
    ba0=d[['bax','bay','baz']].iloc[0].to_numpy();bg0=d[['bgx','bgy','bgz']].iloc[0].to_numpy()
    # Variables are scaled local pose corrections, normalized motion states,
    # constant gyro bias, and constant accelerometer bias.
    ps=.01;bs=.05;gs=.001;N=21*n+6;bgidx=21*n;baidx=bgidx+3
    x=np.zeros(N);x[baidx:]=ba0/bs;x[bgidx:baidx]=bg0/gs
    def decode(x):
        R=q0*Rot.from_rotvec(x[:3*n].reshape(n,3)*ps)
        z=x[3*n:21*n].reshape(n,18);bg=x[bgidx:baidx]*gs;ba=x[baidx:]*bs
        return R,z,bg,ba
    def residual(x):
        R,z,bg,ba=decode(x)
        inc=Rot.identity(n-1)
        for j in range(m):inc=inc*Rot.from_rotvec((rates[:,j]-bg)*hgyro)
        rg=(inc.inv()*(R[:-1].inv()*R[1:])).as_rotvec()/sdgyro
        ra=(R.inv().apply(z[:,[2,8,14]]*scale[[2,8,14]]-g)+ba-acc)/model['acc_sigma']
        rm=(R.inv().apply(np.tile(field,(n,1)))-mag)/model['mag_sigma']
        rp=(z[1:]-z[:-1]@F.T)@Lq.T
        return np.r_[ra.ravel(),rm.ravel(),rg.ravel(),rp.ravel(),Lp@z[0],ba/bias_prior_sd]
    nr=6*n+3*(n-1)+18*(n-1)+21
    pat=sparse.lil_matrix((nr,N),dtype=int)
    for k in range(n):
        pat[3*k:3*k+3,3*k:3*k+3]=1
        for ax in range(3):pat[3*k:3*k+3,3*n+18*k+6*ax+2]=1
        pat[3*k:3*k+3,baidx:]=1;pat[3*n+3*k:3*n+3*k+3,3*k:3*k+3]=1
    off=6*n
    for k in range(n-1):
        pat[off+3*k:off+3*k+3,3*k:3*k+6]=1;pat[off+3*k:off+3*k+3,bgidx:baidx]=1
    off+=3*(n-1)
    for k in range(n-1):pat[off+18*k:off+18*k+18,3*n+18*k:3*n+18*k+36]=1
    pat[-21:-3,3*n:3*n+18]=1
    pat[-3:,baidx:]=1
    structure=pat.tocsr();groups=group_columns(structure)
    def solve(x,steps):
        x=x.copy();history=[];status=0;nfev=0
        for _ in range(steps):
            nfev+=1
            r=residual(x);cost=float(r@r/2)
            J=approx_derivative(residual,x,method='3-point',sparsity=(structure,groups))
            gradient=np.asarray(J.T@r).ravel();optimality=float(np.max(abs(gradient)))
            history.append(dict(cost=cost,gradient=optimality))
            if optimality<1e-5:status=1;break
            normal=(J.T@J).tocsc();ridge=1e-10
            delta=-spsolve(normal+sparse.eye(N,format='csc')*ridge,gradient)
            accepted=False
            for alpha in 2.**-np.arange(15):
                trial=x+alpha*delta;rr=residual(trial);newcost=float(rr@rr/2)
                if newcost<cost:
                    x=trial;accepted=True;break
            if not accepted:
                status=2 if optimality<1e-3 else -1;break
            if np.linalg.norm(alpha*delta)<1e-7*(1+np.linalg.norm(x)):
                status=2;break
        r=residual(x);J=approx_derivative(residual,x,method='3-point',sparsity=(structure,groups))
        return SimpleNamespace(x=x,cost=float(r@r/2),jac=J,optimality=float(np.max(abs(J.T@r))),nfev=nfev,status=status,message='sparse Gauss-Newton with line search; 1e-10 step ridge',history=history)
    first=solve(x,1);fit=solve(first.x,max_nfev)
    def score(x):
        # Only this function reads truth. No truth initializes or enters residuals.
        R,z,bg,ba=decode(x);tq=Rot.from_quat(d[['true_qx','true_qy','true_qz','true_qw']].to_numpy())
        be=ba-d[['true_bax','true_bay','true_baz']].iloc[0].to_numpy()
        return dict(ba_h_mm_s2=float(np.linalg.norm(be[:2])*1000),ba_mm_s2=(be*1000).tolist(),attitude_deg=float(np.sqrt(np.mean((tq.inv()*R).magnitude()**2))*180/np.pi),ba=ba.tolist(),bg=bg.tolist())
    # Schur complement from residual Jacobian: eliminate every motion state,
    # gyro bias and all other attitudes. Scaling: .01 rad and .05 m/s^2.
    J=fit.jac.tocsc()
    # Express initial attitude perturbations about world X/Y/Z. Eliminate
    # world yaw as another nuisance for the five-dimensional tilt/bias test.
    baseR=decode(fit.x)[0][0];S=np.zeros((3,3));eps=1e-4
    for j in range(3):
        xp=fit.x.copy();xm=fit.x.copy();xp[j]+=eps;xm[j]-=eps
        S[:,j]=((decode(xp)[0][0]*baseR.inv()).as_rotvec()-(decode(xm)[0][0]*baseR.inv()).as_rotvec())/(2*eps)
    transform=np.linalg.solve(S,np.eye(3)*ps)
    J=sparse.hstack([J[:,:3]@transform,J[:,3:]],format='csc')
    def information(rows,chosen=None):
        targets=np.r_[np.arange(3),np.arange(baidx,N)] if chosen is None else chosen
        nuisance=np.setdiff1d(np.arange(N),targets)
        A=J[rows,:][:,targets].toarray();B=J[rows,:][:,nuisance];normal=(B.T@B).tocsc()
        # Ridge is a numerical rank regularizer, explicitly reported. Compare
        # two values to expose sensitivity; this is not promoted as physical information.
        vals=[]
        for ridge in (1e-10,1e-12):
            X=spsolve(normal+sparse.eye(normal.shape[0],format='csc')*ridge,B.T@A)
            S=A.T@A-A.T@(B@X);S=(S+S.T)/2
            vals.append(dict(ridge=ridge,eigenvalues=np.linalg.eigvalsh(S).tolist(),diagonal=np.diag(S).tolist()))
        return vals
    out=dict(prefix=str(prefix),nodes=n,dt=h,prior_scale=prior_scale,bias_prior_sd=bias_prior_sd,initial=score(x),one_step=score(first.x),relinearized=score(fit.x),cost=float(fit.cost),optimality=float(fit.optimality),nfev=fit.nfev,status=fit.status,message=fit.message,iterations=fit.history,Q=qdiag,stationary=pdiag,
             motion_information=information(slice(0,nr-3)),full_information=information(slice(None)),measurement_information=information(slice(0,6*n+3*(n-1))),tilt_bias_information=information(slice(None),np.r_[0,1,np.arange(baidx,N)]),
             tilt_bias_motion_information=information(slice(0,nr-3),np.r_[0,1,np.arange(baidx,N)]),
             tilt_bias_measurement_information=information(slice(0,6*n+3*(n-1)),np.r_[0,1,np.arange(baidx,N)]),
             model_stationary_relerr=float(np.linalg.norm(P-np.array(model['stationary18']))/np.linalg.norm(P)),
             target_scales=[ps]*3+[bs]*3,tilt_bias_scales=[ps]*2+[bs]*3)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('prefix',type=Path);p.add_argument('--prior-scale',type=float,default=1);p.add_argument('--max-nfev',type=int,default=30);p.add_argument('--bias-prior-sd',type=float,default=.008*9.80665/np.sqrt(3));p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=run(a.prefix,a.prior_scale,a.max_nfev,a.bias_prior_sd);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
