"""Analytical symmetry reduction for the exact planar MOVING compatibility record.

This proves an invariant coordinate pattern for the ideal real-arithmetic MEKF
mean once the canonical magnetic reference is (B,0,0). It deliberately does
NOT claim an all-time MAGNETIC SERVICE floor: covariance and accepted-service
recurrences remain to be enclosed on the reduced periodic system.
"""
from __future__ import annotations
import numpy as np

def skew(v):
    x,y,z=map(float,v)
    return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def ry(a):
    c,s=np.cos(a),np.sin(a)
    return np.array([[c,0.,s],[0.,1.,0.],[-s,0.,c]])

def planar_record(psi, psi_dot, ax, g=9.80665, cos_theta=39999/40001, B=75.):
    # Delivered record from the exact pair, in the common body frame.
    R=ry(-psi)
    acc=R@np.array([ax,0.,-g*cos_theta])
    gyro=np.array([0.,psi_dot,0.])
    mag=R@np.array([B,0.,0.])
    return acc,gyro,mag

def mean_parity_residuals(psi, psi_dot, ax, aw_xz, ba_xz, g=9.80665,
                          cos_theta=39999/40001, B=75.):
    """Check the coordinate parities used by the invariant-manifold proof.

    Nominal world->body attitude is Ry(-psi); BG_y may evolve, BA/AW/LIN live
    in xz while their y components vanish. The returned identities show that
    accel/mag innovations have zero y component and their attitude Jacobians
    map a pure-y attitude correction back to xz innovations.
    """
    R=ry(-psi)
    acc,gyro,mag=planar_record(psi,psi_dot,ax,g,cos_theta,B)
    aw=np.array([aw_xz[0],0.,aw_xz[1]])
    ba=np.array([ba_xz[0],0.,ba_xz[1]])
    pred_acc=R@(aw-np.array([0.,0.,g]))+ba
    pred_mag=R@np.array([B,0.,0.])
    ra=acc-pred_acc; rm=mag-pred_mag
    Ja=-skew(R@(aw-np.array([0.,0.,g])))
    Jm=-skew(pred_mag)
    ey=np.array([0.,1.,0.])
    return {
      "acc_residual_y":float(ra[1]),"mag_residual_y":float(rm[1]),
      "acc_Jatt_ey_y":float((Ja@ey)[1]),"mag_Jatt_ey_y":float((Jm@ey)[1]),
      "gyro_x":float(gyro[0]),"gyro_z":float(gyro[2]),
      "mag_exact_residual_norm":float(np.linalg.norm(rm)),
    }

def certificate():
    phases=np.linspace(-0.02,0.02,17)
    worst=0.
    for p in phases:
        q=mean_parity_residuals(p,0.003,0.001,[0.002,-0.001],[0.0004,-0.0003])
        worst=max(worst,max(abs(q[k]) for k in ("acc_residual_y","mag_residual_y",
                                                "acc_Jatt_ey_y","mag_Jatt_ey_y",
                                                "gyro_x","gyro_z","mag_exact_residual_norm")))
    return {
      "qualification":"OU3_PLANAR_MOVING_SYMMETRY_V1",
      "result_type":"analytical coordinate-parity identities; sampled arithmetic regression only",
      "mean_invariant_pattern":{
        "attitude":"world-to-body Ry(-psi) on the canonical planar branch",
        "gyro_and_BG":"only y component may be nonzero",
        "delivered_accel_mag":"only x,z components",
        "BA_AW_LIN":"x,z plane; y mean components remain zero under parity-preserving gains",
        "mag_reference":"canonical (B,0,0) after exact symmetric acquisition/refinement",
      },
      "important_limitation":"full covariance is not diagonal/scalar under pitch propagation; magnetic service must retain the reduced attitude/BG covariance and actual S_m recurrence",
      "continuous_hard_iron":"zero-offset planar record is symmetry-compatible, but all-time estimator eligibility/float32 transfer remains a separate obligation",
      "sampled_identity_worst_abs":worst,
      "all_time_magnetic_service_verified":False,
      "theorem_closed":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
