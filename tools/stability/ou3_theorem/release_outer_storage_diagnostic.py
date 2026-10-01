"""Literal first-A21-release outer-storage diagnostic; finite evidence only."""
from pathlib import Path
import json, math, subprocess, tempfile
import numpy as np

REPO=Path(__file__).resolve().parents[3]

def driver_source():
    s=(REPO/"tools/stability/ag_readout_source.cpp").read_text()
    s=s.replace("for (int k=1; k<=45064; ++k) {","for (int k=1; k<=40000; ++k) {")
    s=s.replace("if (k==45001) {","if (k==39950) {")
    s=s.replace("const float roll = wave ? static_cast<float>(.02*std::sin(.5*t)) : 0.0f;","const float roll=0.0f;")
    s=s.replace("const float rate = wave ? static_cast<float>(.01*std::cos(.5*t)) : 0.0f;","const float rate=0.0f;")
    s=s.replace("const float az = wave ? static_cast<float>(-.144*std::sin(.6*t)) : 0.0f;","const float az = wave ? static_cast<float>(6*std::sin(2*t)) : 0.0f;")
    s=s.replace("rwb*Eigen::Vector3f(0,0,az-g_std)","rwb*Eigen::Vector3f(az,0,az-g_std)")
    s=s.replace("rwb*Eigen::Vector3f(60,0,30)","rwb*Eigen::Vector3f(45,0,45)")
    s=s.replace(" || applied!=8","")
    s=s.replace("static_cast<double>(k)*.005","static_cast<double>(k)*static_cast<double>(.005f)")
    s=s.replace("int live=-1, refined=-1, active=-1, applied=0;",
                'int live=-1, refined=-1, active=-1, applied=0; std::string release_cov,release_state,release_quat;')
    old="if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;"
    new='''if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) {
            active=k; const auto& rm=filter.raw().mekf();
            release_cov=matrix_json(rm.Pext); release_state=matrix_json(rm.xext);
            release_quat=matrix_json(rm.qref.coeffs());
        }'''
    if old not in s: raise ValueError("active-release anchor changed")
    s=s.replace(old,new)
    anchor='<< ",\\\"root_covariance\\\":" << root'
    if anchor not in s: raise ValueError("output anchor changed")
    s=s.replace(anchor,'<< ",\\\"release_covariance\\\":" << release_cov'
                       ' << ",\\\"release_state\\\":" << release_state'
                       ' << ",\\\"release_quaternion\\\":" << release_quat '
                       +anchor)
    return s

def rotvec(qcoeff):
    q=np.asarray(qcoeff,dtype=float).reshape(4); xyz=q[:3]; w=float(q[3])
    n=float(np.linalg.norm(xyz))
    if n==0: return np.zeros(3)
    angle=2*math.atan2(n,w)
    if angle>math.pi: angle-=2*math.pi
    return xyz*(angle/n)

def audit(native):
    P=np.asarray(native["release_covariance"],dtype=float)
    theta=rotvec(native["release_quaternion"])
    Ptt=P[:3,:3]
    eig=np.linalg.eigvalsh(Ptt)
    if eig[0]<=0: raise ArithmeticError("attitude marginal is not SPD")
    Vtheta=float(theta@np.linalg.solve(Ptt,theta))
    return {
      "active_step":int(native["active_step"]),
      "release_time_s":int(native["active_step"])*float(np.float32(.005)),
      "release_tilt_rad":float(np.linalg.norm(theta)),
      "release_tilt_deg":float(np.linalg.norm(theta))*180/math.pi,
      "release_attitude_covariance":Ptt.tolist(),
      "release_attitude_cov_eigenvalues":eig.tolist(),
      "ba_eliminated_outer_storage_lower_from_attitude":Vtheta,
      "ba_eliminated_outer_sqrt_storage_lower_from_attitude":math.sqrt(Vtheta),
      "target_radius":.15,"target_squared_radius":.0225,
      "attitude_alone_excludes_entry":Vtheta>.0225,
      "qualification":"FINITE_CARRIED_RELEASE_DIAGNOSTIC_ONLY",
      "source_uniform_verified":False}

def run(eigen=Path("/usr/include/eigen3")):
    with tempfile.TemporaryDirectory(prefix="ou3-release-outer-") as d:
        src=Path(d)/"release.cpp"; exe=Path(d)/"release"; src.write_text(driver_source())
        subprocess.run(["g++","-O2","-std=c++20","-I"+str(REPO/"src"),"-isystem",str(eigen),str(src),"-o",str(exe)],check=True)
        native=json.loads(subprocess.check_output([str(exe),"wave"],text=True))
    return {"audit":audit(native),"native_release":{k:native[k] for k in
      ("active_step","release_covariance","release_state","release_quaternion")}}

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("--eigen",type=Path,default=Path("/usr/include/eigen3")); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); a.output.write_text(json.dumps(run(a.eigen),indent=2,sort_keys=True)+"\n")
