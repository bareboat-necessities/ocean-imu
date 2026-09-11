"""Conditional shipping graph, not a reachable source cover or stability gate.

Operands are functions of the SAME pre-operation state. H is derived from the
measurement model, never from P^-1 N. Float execution defects remain explicit;
these real-arithmetic identities alone do not enclose deployment roundoff.
"""
from __future__ import annotations
import math
import numpy as np


def skew(v):
    x,y,z=np.asarray(v,dtype=float)
    return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])


def quaternion_matrix(q):
    """Eigen's quaternion-to-matrix polynomial, without silently normalizing q."""
    w,x,y,z=np.asarray(q,dtype=float)
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])


def measurement_graph(kind,P,q,x,measured,noise,active,mag_reference,gravity):
    """Real-arithmetic N/S/r for zero lever arm and temperature-centered bias.

    This is a conditional identity for all operands in that branch, not a
    selection of independent N/S boxes. Shipping builds S from the displayed
    upper cross blocks and their transposes; retaining that evaluation matters
    when the recorded binary32 P has a small asymmetry. H_residual is full even
    when H18 masks the gain's bias columns and rows.
    """
    P=np.asarray(P,dtype=float).reshape(21,21)
    x=np.asarray(x,dtype=float).reshape(21)
    noise=np.asarray(noise,dtype=float).reshape(3,3)
    measured=np.asarray(measured,dtype=float).reshape(3)
    rot=quaternion_matrix(q)
    H=np.zeros((3,21))
    if kind==10:
        cog=rot@(x[15:18]-np.array([0.,0.,gravity]))
        J=-skew(cog); A=rot
        H[:,:3]=J;H[:,15:18]=A;H[:,18:21]=np.eye(3)
        r=measured-cog-x[18:21]
        S=(noise+J@P[:3,:3]@J.T+J@P[:3,15:18]@A.T
           +A@P[:3,15:18].T@J.T+A@P[15:18,15:18]@A.T
           +J@P[:3,18:21]+P[:3,18:21].T@J.T
           +A@P[15:18,18:21]+P[15:18,18:21].T@A.T+P[18:21,18:21])
    elif kind==11:
        predicted=rot@np.asarray(mag_reference,dtype=float)
        H[:,:3]=-skew(predicted)
        r=measured-predicted
        S=noise+H[:,:3]@P[:3,:3]@H[:,:3].T
    elif kind==12:
        H[:,12:15]=np.eye(3);r=-x[12:15]
        S=P[12:15,12:15]+noise
    else:
        raise ValueError('not a supported measurement event')
    Hg=H.copy()
    if kind==10 and not active:Hg[:,18:21]=0
    N=P@Hg.T
    if not active:N[18:21]=0
    return {'N':N,'S':S,'r':r,'H_residual':H,'H_gain':Hg}


def binding_defects(row):
    params=row['parameters'].astype(float)
    if params[3]!=0:
        raise ValueError('lever extension has no analytic binding here')
    graph=measurement_graph(int(row['kind']),row['P'],row['q'],row['x'],
        row['measured'],row['R'],bool(row['active']),row['mag_reference'],params[0])
    result={}
    for name in ('N','S','r'):
        actual=row[name].astype(float).reshape(graph[name].shape)
        result[name+'_absolute_defect']=float(np.max(np.abs(actual-graph[name])))
        result[name+'_scaled_defect']=result[name+'_absolute_defect']/max(1.,float(np.max(np.abs(actual))))
    # Current probe has zero de-heel. Do not assume that if the default changes.
    if params[4]!=0:raise ValueError('nonzero de-heel needs a physical input binding')
    kind=int(row['kind'])
    raw=(row['source'][15:18].astype(float) if kind==10
         else row['source'][21:24].astype(float)-row['hard_iron'].astype(float))
    result['physical_input_max_abs_defect']=(float(np.max(np.abs(row['measured']-raw))) if kind!=12 else 0.)
    result['held_bias_N_max_abs']=float(np.max(np.abs(row['N'].reshape(21,3)[18:21]))) if not row['active'] else 0.
    return result


def same_driver_bias_graph(beta_before,beta_after,bhat_before,bhat_after,phi_true,phi_hat):
    """One w, not two independent driver slots; mismatch stays in the joint map."""
    beta0,beta1,bh0,bh1=map(lambda v:np.asarray(v,dtype=float),
                           (beta_before,beta_after,bhat_before,bhat_after))
    w=beta1-phi_true*beta0
    e0,e1=beta0-bh0,beta1-bh1
    expected=phi_hat*e0+(phi_true-phi_hat)*beta0+w
    return w,e1-expected


def binding_report(rows,family,dt):
    """Observed residuals and event census. No maximum is a universal bound."""
    maxima={};bias_max=0.;driver_max=0.;bias_count=0;solve_count=0
    mean_max=joseph_max=reset_max=projection_max=0.
    last_sample=None;pre_prediction=None;last_solve=None;last_joseph=None;pre_projection=None
    for row in rows:
        kind=int(row['kind'])
        if kind in (10,11,12):
            for name,value in binding_defects(row).items():maxima[name]=max(maxima.get(name,0.),value)
            last_solve=row;solve_count+=1
        elif kind in (50,51,52):
            if last_solve is None or int(last_solve['kind'])!=kind-40:raise ValueError('mean event ancestry mismatch')
            before=last_solve['x'].astype(float)
            expected=before+last_solve['K'].reshape(21,3).astype(float)@last_solve['r'].astype(float)
            mean_max=max(mean_max,float(np.max(np.abs(row['x']-expected))))
        elif kind in (60,61,62):
            if last_solve is None or int(last_solve['kind'])!=kind-50:raise ValueError('Joseph event ancestry mismatch')
            p=last_solve['P'].reshape(21,21).astype(float)
            n=last_solve['N'].reshape(21,3).astype(float);k=last_solve['K'].reshape(21,3).astype(float)
            s=last_solve['S'].reshape(3,3).astype(float)
            raw=p-k@n.T-n@k.T+k@s@k.T
            # Literal upper-triangle write and mirror, not an independent symmetric S.
            expected=np.triu(raw)+np.triu(raw,1).T
            joseph_max=max(joseph_max,float(np.max(np.abs(row['P'].reshape(21,21)-expected))))
            last_joseph=row
        elif kind==70:
            if last_joseph is None:raise ValueError('reset without Joseph ancestor')
            g=np.eye(21);g[:3,:3]+=.5*skew(last_joseph['x'][:3])
            expected=g@last_joseph['P'].reshape(21,21).astype(float)@g.T
            reset_max=max(reset_max,float(np.max(np.abs(row['P'].reshape(21,21)-expected))))
            pre_projection=row
        elif kind==80:
            if pre_projection is None:raise ValueError('projection without reset ancestor')
            b=pre_projection['x'][18:21].astype(float);limit=float(row['parameters'][2])
            expected=b*min(1.,limit/np.linalg.norm(b)) if limit>0 and np.linalg.norm(b)>0 else b
            projection_max=max(projection_max,float(np.max(np.abs(row['x'][18:21]-expected))))
        elif kind==40:last_sample=row
        elif kind==1:pre_prediction=row
        elif kind==2:
            if pre_prediction is None or last_sample is None:raise ValueError('prediction source ancestry missing')
            if int(row['sample'])!=int(last_sample['sample'])+1:raise ValueError('nonconsecutive physical source')
            phi_true=1. if family==2 else math.exp(-dt/(900. if family==0 else 1200.))
            phi_hat=(math.exp(-dt/float(row['parameters'][1])) if pre_prediction['active'] else 1.)
            w,defect=same_driver_bias_graph(last_sample['source'][18:21],row['source'][18:21],
                pre_prediction['x'][18:21],row['x'][18:21],phi_true,phi_hat)
            bias_max=max(bias_max,float(np.max(np.abs(defect))))
            driver_max=max(driver_max,float(np.max(np.abs(w))));bias_count+=1
    # Diagnostic tolerances catch broken attachment, not theorem thresholds.
    defects_ok=(all(maxima.get(k+'_scaled_defect',math.inf)<2e-5 for k in ('N','S','r'))
                and maxima.get('physical_input_max_abs_defect',math.inf)<1e-5)
    return {'qualification':'OBSERVED_OPERAND_BINDING_DEFECTS_NOT_UNIFORM_ENCLOSURE',
        'solves':solve_count,'measurement_graph':maxima,
        'mean_update_max_abs_defect':mean_max,'Joseph_max_abs_defect':joseph_max,
        'reset_covariance_max_abs_defect':reset_max,'projection_max_abs_defect':projection_max,
        'same_driver_bias_predictions':bias_count,'same_driver_bias_graph_max_abs_defect':bias_max,
        'observed_bias_driver_component_max_abs':driver_max,
        'measurement_binding_regression_pass':defects_ok,
        'source_uniform_attachment_certified':False,'roundoff_bounds_certified':False}
