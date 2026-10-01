"""Literal release snapshot and BA-eliminated outer-entry diagnostic.

Finite carried diagnostic only. It snapshots the unchanged shipping execution at
first A21 BA release and computes a rigorous *component lower bound* on the
BA-eliminated outer storage from attitude alone. It is not source-uniform entry.
"""
from pathlib import Path
import json, math, subprocess, tempfile
import numpy as np
from .construction_history_diagnostic import REPO\n\ndef driver_source():\n    return (REPO/"tools/stability/ag_readout_source.cpp").read_text()

def release_driver_source():
    s=driver_source()\n    s=s.replace("for (int k=1; k<=45064; ++k) {", "for (int k=1; k<=40000; ++k) {")\n    s=s.replace("if (k==45001) {", "if (k==39950) {")
    # This probe needs only the first A21 release (~step 36008), not the
    # construction diagnostic's 600-s tail. Keep a short post-release scoring
    # tail so the inherited JSON remains finite.
    s=s.replace("k<=120000", "k<=37000")
    s=s.replace("k>80000", "k>36000")
    s=s.replace("std::string root;", 'std::string root="[]";')
    old="int live=-1, refined=-1, active=-1, applied=0;"
    # driver_source inherits this declaration from ag_readout_source.cpp.
    if old not in s: raise ValueError("release declaration anchor changed")
    s=s.replace(old, old+"\n    std::string release_cov, release_state, release_quat;")
    old_active="if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;"
    new_active='''if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) {
                active=k;
                const auto& rm=filter.raw().mekf();
                release_cov=matrix_json(rm.Pext);
                release_state=matrix_json(rm.xext);
                release_quat=matrix_json(rm.qref.coeffs());
            }'''
    if old_active not in s: raise ValueError("release event anchor changed")
    s=s.replace(old_active,new_active)
    anchor='<< ",\\\"live_step\\\":" << live'
    if anchor not in s: raise ValueError("JSON output anchor changed")
    s=s.replace(anchor, '<< ",\\\"release_covariance\\\":" << release_cov'
                      ' << ",\\\"release_state\\\":" << release_state'
                      ' << ",\\\"release_quaternion\\\":" << release_quat '
                      +anchor)
    return s

def _quat_rotvec(coeffs):
    q=np.asarray(coeffs,dtype=float).reshape(4)
    xyz=q[:3]; w=q[3]
    n=float(np.linalg.norm(xyz))
    if n==0: return np.zeros(3)
    angle=2.0*math.atan2(n,w)
    if angle>math.pi: angle-=2*math.pi
    return xyz/n*angle

def audit(native):
    P=np.asarray(native["release_covariance"],dtype=float)
    if P.shape!=(21,21): raise ValueError("21x21 release covariance required")
    # Eliminating BA from the covariance-metric quadratic means minimizing over
    # e_ba. At release P_ob=0 by the held-BA decoupling, so the outer marginal
    # is also the release Schur block. To identify an obstruction without
    # reconstructing physical v/p/S truth, minimize the outer quadratic further
    # over all non-attitude outer errors. This leaves theta' P_tt^-1 theta.
    theta=_quat_rotvec(native["release_quaternion"])
    Ptt=P[:3,:3]
    eig=np.linalg.eigvalsh(Ptt)
    if eig[0] <= 0: raise ArithmeticError("release attitude marginal not SPD")
    vtheta=float(theta @ np.linalg.solve(Ptt,theta))
    rtheta=math.sqrt(vtheta)
    return {
      "active_step":native["active_step"],
      "release_time_s":native["active_step"]*float(np.float32(.005)),
      "release_tilt_rad":float(np.linalg.norm(theta)),
      "release_tilt_deg":float(np.linalg.norm(theta))*180/math.pi,
      "release_attitude_covariance":Ptt.tolist(),
      "release_attitude_cov_eigenvalues":eig.tolist(),
      "ba_eliminated_outer_storage_lower_from_attitude":vtheta,
      "ba_eliminated_outer_sqrt_storage_lower_from_attitude":rtheta,
      "target_radius":.15,
      "target_squared_radius":.0225,
      "attitude_alone_excludes_entry":vtheta>.0225,
      "qualification":"FINITE_CARRIED_RELEASE_DIAGNOSTIC_ONLY",
      "source_uniform_verified":False
    }

def run(eigen=Path("/usr/include/eigen3")):
    source=release_driver_source()
    with tempfile.TemporaryDirectory(prefix="ou3-release-outer-") as directory:
        driver=Path(directory)/"release.cpp"; binary=Path(directory)/"release"
        driver.write_text(source)
        subprocess.run(["g++","-O2","-std=c++20","-I"+str(REPO/"src"),
                        "-isystem",str(eigen),str(driver),"-o",str(binary)],check=True)
        native=json.loads(subprocess.check_output([str(binary),"wave"],text=True))
    return {"audit":audit(native),"native_release":{
      k:native[k] for k in ("active_step","release_covariance","release_state","release_quaternion")
    }}

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("--eigen",type=Path,default=Path("/usr/include/eigen3"))
    p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    a.output.write_text(json.dumps(run(a.eigen),indent=2,sort_keys=True)+"\n")
