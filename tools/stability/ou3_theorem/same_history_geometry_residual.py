# ruff: noqa: F401, F811
"""Same-history geometry/residual propagation for one release boundary.

The estimator error coordinates are primary.  We do not reconstruct an
independent absolute attitude box: given the physical B->W rotation and the
same-boundary left attitude error theta, shipping W->B geometry follows by
composition.  Residual and reset injection are then generated from the same
state/covariance boundary.
"""
from __future__ import annotations
from dataclasses import dataclass
import math,numpy as np
from .accel_geometry_cell import rodrigues
G=9.80665

def skew(v):
 x,y,z=np.asarray(v,float);return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def physical_rotation_from_gravity_heading(g_body,heading_axis_body=None):
 """Construct B->W from gravity direction plus a supplied horizontal heading.

 Gravity alone leaves yaw free, so source-uniform construction must supply a
 heading/reference direction or use a world-frame attitude-free row theorem.
 """
 z=np.asarray(g_body,float);nz=np.linalg.norm(z)
 if nz<=0:raise ArithmeticError("gravity direction not separated from zero")
 z=z/nz
 if heading_axis_body is None:
  raise ArithmeticError("gravity direction does not determine yaw")
 x=np.asarray(heading_axis_body,float);x=x-z*(z@x);nx=np.linalg.norm(x)
 if nx<=0:raise ArithmeticError("heading parallel to gravity")
 x=x/nx;y=z@np.zeros(3) if False else np.cross(z,x)
 # columns are world axes represented in body; transpose maps body to world.
 R_wb=np.column_stack((x,y,z))
 return R_wb.T

def estimator_R_wb_from_error(R_bw_phys,theta_error):
 """Shipping qref is W->B. true W->B=Exp(theta)*estimated W->B."""
 R_bw=np.asarray(R_bw_phys,float)
 R_true_wb=R_bw.T
 return rodrigues(-np.asarray(theta_error,float))@R_true_wb

def acc_geometry_and_residual(*,R_bw_phys,theta_error,physical_a,
                              accel_error,aw_hat,ba_hat,lever=None,g=G):
 Rhat_wb=estimator_R_wb_from_error(R_bw_phys,theta_error)
 aw=np.asarray(aw_hat,float);a=np.asarray(physical_a,float)
 lever=np.zeros(3) if lever is None else np.asarray(lever,float)
 ba=np.asarray(ba_hat,float);ea=np.asarray(accel_error,float)
 gw=np.array([0.,0.,g])
 # Physical ideal specific force in body coordinates, plus delivered error.
 f_meas=R_bw_phys.T@(a-gw)+ea
 f_cog=Rhat_wb@(aw-gw)
 f_pred=f_cog+lever+ba
 r=f_meas-f_pred
 Jatt=-skew(f_cog);Jaw=Rhat_wb
 return {"R_wb":Rhat_wb,"f_meas":f_meas,"f_cog":f_cog,
         "residual":r,"Jatt":Jatt,"Jaw":Jaw}

def local_measurement_defect(residual,H,error):
 return np.asarray(residual,float)-np.asarray(H,float)@np.asarray(error,float)

def s_geometry_and_residual(S_phys,S_hat):
 r=-np.asarray(S_hat,float)
 e=np.asarray(S_phys,float)-np.asarray(S_hat,float)
 return {"residual":r,"error_S":e,"physical_source":np.asarray(S_phys,float),
         "identity_error":float(np.linalg.norm(r-(e-np.asarray(S_phys,float))))}

def correction_injection(K,residual):
 K=np.asarray(K,float);r=np.asarray(residual,float)
 d=K@r
 return {"state_increment":d,"dtheta":d[:3].copy(),
         "same_gain_residual_product":True}

def boundary(*,R_bw_phys,theta_error,physical_a,accel_error,aw_hat,ba_hat,
             H_acc,error_state,K_acc,S_phys,S_hat,lever=None):
 a=acc_geometry_and_residual(R_bw_phys=R_bw_phys,theta_error=theta_error,
    physical_a=physical_a,accel_error=accel_error,aw_hat=aw_hat,
    ba_hat=ba_hat,lever=lever)
 a["local_defect"]=local_measurement_defect(a["residual"],H_acc,error_state)
 a["injection"]=correction_injection(K_acc,a["residual"])
 return {"acc":a,"S":s_geometry_and_residual(S_phys,S_hat),
         "same_history_boundary":True}

def world_acc_row(*,aw_hat,g=G):
 """Attitude-free world-frame accelerometer attitude geometry.

 By Lemma W, body row -[R f]x factors orthogonally to [f_world]x with
 f_world=aw_hat-g e_z.  This is the preferred source-uniform representation
 when HistoryCell does not carry absolute yaw.
 """
 fw=np.asarray(aw_hat,float)-np.array([0.,0.,g])
 return {"force_world":fw,"Jatt_world":skew(fw),
         "absolute_attitude_required":False,
         "orthogonal_body_factor_suppressed":True}

def world_acc_residual_source(*,physical_a,aw_hat,world_sensor_error,
                              ba_world=None,lever_world=None):
 """World-coordinate residual source before the orthogonal body factor.

 This retains the signed physical-minus-nominal acceleration source.  BA and
 lever terms must already be transported by the same estimator boundary; they
 are not independently oriented here.
 """
 z=np.asarray(physical_a,float)-np.asarray(aw_hat,float)+np.asarray(world_sensor_error,float)
 if ba_world is not None:z=z-np.asarray(ba_world,float)
 if lever_world is not None:z=z-np.asarray(lever_world,float)
 return z

def world_boundary(*,physical_a,aw_hat,world_sensor_error,error_state,H_world,
                   K_world,S_phys,S_hat,ba_world=None,lever_world=None):
 geom=world_acc_row(aw_hat=aw_hat)
 r=world_acc_residual_source(physical_a=physical_a,aw_hat=aw_hat,
       world_sensor_error=world_sensor_error,ba_world=ba_world,lever_world=lever_world)
 return {"acc":{"force_world":geom["force_world"],"Jatt_world":geom["Jatt_world"],
                "residual_world":r,
                "local_defect_world":local_measurement_defect(r,H_world,error_state),
                "injection":correction_injection(K_world,r)},
         "S":s_geometry_and_residual(S_phys,S_hat),
         "same_history_boundary":True,"absolute_yaw_eliminated":True}
