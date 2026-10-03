"""Schur/conditional-information identity for the literal accelerometer update.

For a joint Gaussian prior on target x=(attitude,BG) and nuisance n=(AW,BA,...),
the literal update y=A x+B n+noise changes the target information by the
conditional measurement information after eliminating n. This module proves the
matrix identity numerically/algebraically and records the comparison direction
needed by the planar service oracle.
"""
from __future__ import annotations
import numpy as np

def target_posterior(Pxx,Pxn,Pnn,A,B,R):
    Pxx=np.asarray(Pxx,float);Pxn=np.asarray(Pxn,float);Pnn=np.asarray(Pnn,float)
    A=np.asarray(A,float);B=np.asarray(B,float);R=np.asarray(R,float)
    C=np.concatenate((A,B),axis=1)
    P=np.block([[Pxx,Pxn],[Pxn.T,Pnn]])
    S=C@P@C.T+R
    K=P@C.T@np.linalg.inv(S)
    Q=P-K@S@K.T
    return Q[:Pxx.shape[0],:Pxx.shape[0]]

def conditional_information(Pxx,Pxn,Pnn,A,B,R):
    # n|x has covariance Pnn-Pnx Pxx^-1 Pxn and mean slope Pnx Pxx^-1 x.
    ix=np.linalg.inv(Pxx)
    L=Pxn.T@ix
    N=Pnn-Pxn.T@ix@Pxn
    Ae=A+B@L
    Re=R+B@N@B.T
    return Ae.T@np.linalg.inv(Re)@Ae

def identity_defect(Pxx,Pxn,Pnn,A,B,R):
    post=target_posterior(Pxx,Pxn,Pnn,A,B,R)
    info=np.linalg.inv(post)-np.linalg.inv(Pxx)
    cond=conditional_information(Pxx,Pxn,Pnn,A,B,R)
    return float(np.linalg.norm(info-cond)),info,cond

def oracle_dominates_literal_if(info_oracle,info_literal,tol=1e-12):
    # "More informative oracle" means I_oracle >= I_literal. This is the
    # direction relevant to a worst-case magnetic-service comparison: it can
    # only shrink target covariance more before the next mag sample.
    d=np.asarray(info_oracle,float)-np.asarray(info_literal,float)
    return float(np.linalg.eigvalsh((d+d.T)/2).min())>=-tol

def certificate():
    Pxx=np.array([[.002,.0001],[.0001,.001]])
    Pxn=np.array([[.0002,-.0001],[.00005,.00008]])
    Pnn=np.array([[.003,.0002],[.0002,.002]])
    A=np.array([[2.,.1],[0.,1.]])
    B=np.array([[1.,0.],[0.,.5]])
    R=.04*np.eye(2)
    defect,_,_=identity_defect(Pxx,Pxn,Pnn,A,B,R)
    return {"qualification":"OU3_ACCEL_CONDITIONAL_INFORMATION_V1",
            "schur_identity_defect":defect,
            "identity":"J_x^+-J_x^- = A_eff^T R_eff^-1 A_eff after exact Gaussian nuisance elimination",
            "comparison_needed":"prove I_acc,literal <= I_acc,oracle samplewise on the same planar cell; then establish that this ordering is sufficient for the magnetic-service lower comparison through prediction/reset/mag chronology",
            "important_warning":"samplewise accelerometer information ordering alone is not yet a proof of terminal magnetic Gramian ordering because homogeneous probes and covariance are both updated",
            "shipping_service_proved":False}

if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
