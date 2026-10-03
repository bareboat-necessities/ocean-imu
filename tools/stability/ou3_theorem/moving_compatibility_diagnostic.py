"""Finite literal shipping replay of the exact moving compatibility family.

Read-only instrumentation exports actual F, H, K, final innovation S and reset G.
No estimator/tuner/frontend/scheduler/source state is restarted at service roots.
Only homogeneous diagnostic columns are restarted. This samples disjoint 1-s
windows; it does NOT certify every placed window or any all-time continuation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[3]
HEADER = REPO / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
PROBE = Path(__file__).with_name("moving_compatibility_probe.cpp")

def instrument(text):
    s = text
    def once(a,b):
     nonlocal s
     assert s.count(a)==1,(a,s.count(a));s=s.replace(a,b)
    once('    apply_pending_aw_covariance_inflation_();', '''    if (moving_recording) {
            Eigen::Matrix<T,3,3> Qba_diag = Eigen::Matrix<T,3,3>::Zero();
            if constexpr (with_accel_bias) {
                const T tau_b_dbg = std::max(T(1e-3), tau_bacc_);
                const T phi_b_dbg = acc_bias_updates_enabled_ ? std::exp(-Ts/tau_b_dbg) : T(1);
                if (acc_bias_updates_enabled_) {
                    const T qd_dbg = -T(0.5)*tau_b_dbg*std::expm1(-T(2)*Ts/tau_b_dbg);
                    Qba_diag = Q_bacc_*qd_dbg;
                }
                moving_prediction(F_AA.data(),F_LL.data(),Q_AA.data(),Q_LL.data(),phi_b_dbg,Qba_diag.data(),Pext.data());
            } else {
                moving_prediction(F_AA.data(),F_LL.data(),Q_AA.data(),Q_LL.data(),T(1),Qba_diag.data(),Pext.data());
            }
        }
        apply_pending_aw_covariance_inflation_();''')
    for start,stop,kind,hcode in [
     ('::measurement_update_acc_only(','::measurement_update_mag_only(','acc','H.template block<3,3>(0,0)=J_att; H.template block<3,3>(0,3)=J_bg; H.template block<3,3>(0,15)=J_aw; H.template block<3,3>(0,18).setIdentity();'),
     ('::measurement_update_mag_only(','::accelerometer_measurement_func(','mag','H.template block<3,3>(0,0)=J_att;'),
     ('::applyIntegralZeroPseudoMeas()','::measurement_update_position_pseudo(','S','H.template block<3,3>(0,12).setIdentity();')]:
     a,b=s.index(start),s.index(stop);part=s[a:b];assert part.count('    xext.noalias() += K * r;')==1
     part=part.replace('    xext.noalias() += K * r;',f'''    if (moving_recording) {{
            Eigen::Matrix<T,3,NX> H=Eigen::Matrix<T,3,NX>::Zero(); {hcode}
            moving_correction("{kind}",H.data(),S_mat.data(),K.data());
        }}
        xext.noalias() += K * r;''')
     s=s[:a]+part+s[b:]
    once('    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);', '''    if (moving_recording) moving_reset(dtheta_injected.data());
        ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);''')
    return s

def run_native(eigen, duration=240.0):
    eigen = Path(eigen)
    if not (eigen / "Eigen/Dense").is_file():
        raise ValueError("Eigen include directory must contain Eigen/Dense")
    if not 220 <= duration <= 3600:
        raise ValueError("duration must lie in [220,3600] seconds")
    tapped = instrument(HEADER.read_text())
    with tempfile.TemporaryDirectory(prefix="ou3-moving-compatibility-") as directory:
        work = Path(directory)
        path = work / "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
        path.parent.mkdir(parents=True)
        path.write_text(tapped)
        binary = work / "probe"
        command = ["g++", "-std=c++20", "-O1", "-DEIGEN_UNROLLING_LIMIT=0", "-DEIGEN_NON_ARDUINO",
                   "-I"+str(work), "-I"+str(REPO/"src"), "-I"+str(eigen),
                   str(PROBE), "-o", str(binary)]
        p = subprocess.run(command, capture_output=True, text=True, timeout=90)
        if p.returncode:
            raise RuntimeError("probe compile failed:\\n" + p.stdout + "\\n" + p.stderr)
        result = subprocess.run([str(binary), str(duration)], check=True,
                                capture_output=True, text=True, timeout=240)
        native = json.loads(result.stdout)
        def check_finite(value):
            if isinstance(value, dict):
                return all(check_finite(v) for v in value.values())
            if isinstance(value, list):
                return all(check_finite(v) for v in value)
            return not isinstance(value, (int, float)) or math.isfinite(value)
        if not check_finite(native):
            raise ValueError("native diagnostic emitted a nonfinite quantity")
        if native["service_windows"] <= 0 or native["service_mags_min"] <= 0:
            raise ValueError("native diagnostic contains no applied magnetic service windows")
    return {
        "qualification": "OU3_MOVING_COMPATIBILITY_CARRIED_V1",
        "result_type": "finite carried float32 diagnostic; not an interval theorem",
        "native": native,
        "probe_sha256": hashlib.sha256(PROBE.read_bytes()).hexdigest(),
        "shipping_header_sha256": hashlib.sha256(HEADER.read_bytes()).hexdigest(),
        "instrumented_header_sha256": hashlib.sha256(tapped.encode()).hexdigest(),
        "noise_profile_source": "moving_compatibility_probe.cpp: sigma_a=0.2 matches current theorem/default wrapper; sigma_g=0.00135 and sigma_m=0.8 match the world-frame fixture; this is not that fixture unchanged",
        "same_as_world_frame_fixture": False,
        "noise_profile": {"sigma_a":0.2,"sigma_g":0.00135,"sigma_m":0.8},
        "accepted_information_used": True,
        "actual_innovation_covariance_used": True,
        "actual_gain_and_reset_chronology_used": True,
        "service_window_coverage": "disjoint one-second sampled-root windows after BA activation plus 17 s; additionally every IMU-sample root over the final 20-s source phase cycle",
        "reduced_covariance_diagnostic": "tail P_theta eigenvalue envelope plus BG and theta-BG block norms; finite carried feasibility only",
        "all_placed_service_windows_verified": False,
        "all_time_magnetic_service_verified": False,
        "nonlinear_or_float32_all_time_totality_verified": False,
        "full_shipping_counterexample_admitted": False,
        "shipping_behavior_changed": False,
        "theorem_closed": False,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eigen",type=Path,default=Path("/usr/include/eigen3"))
    p.add_argument("--duration",type=float,default=240.0)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    text=json.dumps(run_native(args.eigen,args.duration),indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text)
    else:
        print(text,end="")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
