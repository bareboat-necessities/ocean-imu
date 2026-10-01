"""Literal A21-release snapshot and BA-eliminated outer-storage diagnostic."""
from pathlib import Path
import json, math, subprocess, tempfile
import numpy as np
from .construction_history_diagnostic import REPO

def driver_source():
    s=(REPO/"tools/stability/ag_readout_source.cpp").read_text()
    s=s.replace("for (int k=1; k<=45064; ++k) {","for (int k=1; k<=37000; ++k) {")
    s=s.replace("if (k==45001) {","if (k==36950) {")
    s=s.replace("int live=-1, refined=-1, active=-1, applied=0;",
      "int live=-1, refined=-1, active=-1, applied=0; std::string release_cov,release_state,release_quat;")
    old="if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;"
    new="""if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) {
            active=k; const auto& rm=filter.raw().mekf();
            release_cov=matrix_json(rm.Pext); release_state=matrix_json(rm.xext);
            release_quat=matrix_json(rm.qref.coeffs());
        }"""
    if old not in s: raise ValueError("active anchor changed")
    s=s.replace(old,new)
    oldout='<< ",\\\"root_covariance\\\":" << root'
    newout='<< ",\\\"release_covariance\\\":" << release_cov << ",\\\"release_state\\\":" << release_state << ",\\\"release_quaternion\\\":" << release_quat << ",\\\"root_covariance\\\":" << root'
    if oldout not in s: raise ValueError("output anchor changed")
    return s.replace(oldout,newout)

def rotvec(coeffs):
    q=np.asarray(coeffs,float).reshape(4); xyz=q[:3]; w=q[3]; n=np.linalg.norm(xyz)
    if n==0: return np.zeros(3)
    a=2*math.atan2(n,w)
    if a>math.pi: a-=2*math.pi
    return xyz/n*a

def audit(native):
    P=np.asarray(native["release_covariance"],float)
    theta=rotvec(native["release_quaternion"])
    Ptt=P[:3,:3]; eig=np.linalg.eigvalsh(Ptt)
    val=float(theta@np.linalg.solve(Ptt,theta))
    return {"active_step":native["active_step"],
      "release_time_s":native["active_step"]*float(np.float32(.005)),
      "release_tilt_deg":float(np.linalg.norm(theta)*180/math.pi),
      "attitude_storage_lower_after_BA_and_all_other_outer_elimination":val,
      "attitude_sqrt_storage_lower":math.sqrt(val),
      "attitude_cov_eigenvalues":eig.tolist(),"target_radius":.15,
      "target_squared_radius":.0225,"attitude_alone_excludes_entry":val>.0225,
      "qualification":"FINITE_CARRIED_RELEASE_DIAGNOSTIC_ONLY","source_uniform_verified":False}

def run(eigen):
    with tempfile.TemporaryDirectory(prefix="ou3-release-") as d:
        src=Path(d)/"r.cpp"; exe=Path(d)/"r"; src.write_text(driver_source())
        subprocess.run(["g++","-O2","-std=c++20","-I"+str(REPO/"src"),"-isystem",str(eigen),str(src),"-o",str(exe)],check=True)
        native=json.loads(subprocess.check_output([str(exe),"wave"],text=True))
    return {"audit":audit(native),"native_release":{k:native[k] for k in ("active_step","release_covariance","release_state","release_quaternion")}}

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("--eigen",type=Path,default=Path("/usr/include/eigen3")); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); a.output.write_text(json.dumps(run(a.eigen),indent=2,sort_keys=True)+"\n")
