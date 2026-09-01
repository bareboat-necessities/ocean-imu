"""Check physical error arithmetic against the shipping C++ filter."""
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval, matrix_identity
import ou3_interval_ad as AD
import ou3_p4_h18_differential_operations as OPS
import ou3_p5_full_h_prefix_cells as FULL
import ou3_source_reachable_matrix_p3 as P3
import ou3_p4_source_path_reachability as PATH


def state(values=None):
    z = [0.0]*18
    for i, x in (values or {}).items():
        z[i] = x
    return AD.independent_vector([Interval.point(x) for x in z])


def domain():
    return json.loads((ROOT/"tools/ou3_proof_operating_domain.json").read_text())


def s_update():
    z = state({0: 0.4, 12: 1.0})
    P = matrix_identity(18)
    P[0][0] = Interval.point(2.0)
    P[1][12] = P[12][1] = Interval.point(0.01)
    return OPS.accepted_update(
        P, z, FULL._H_S(), matrix_identity(3), OPS.S_residual(z))


class SourceMapParityTests(unittest.TestCase):
    def setUp(self):
        self.saved_n = FULL.N
        FULL.N = 18

    def tearDown(self):
        FULL.N = self.saved_n

    def encloses(self, iv, value, tol=0.0):
        self.assertLessEqual(iv.lo, value+tol)
        self.assertGreaterEqual(iv.hi, value-tol)

    def test_S_error_update_contracts_in_its_Joseph_metric(self):
        z = state({12: 1.0})
        post, out, _ = OPS.accepted_update(
            matrix_identity(18), z, FULL._H_S(), matrix_identity(3), OPS.S_residual(z))
        self.encloses(out[12].val, 0.5)
        self.assertFalse(out[12].val.contains(1.5))
        self.encloses(out[12].der[12], 0.5)
        self.encloses(post[12][12], 0.5)
        self.assertLess(out[12].val.hi**2/post[12][12].lo, 0.500000000001)

    def test_wave_truth_S_is_an_external_input(self):
        z = state({12: 3.0})
        r = OPS.S_residual(z, truth_S=[Interval.point(x) for x in (5,0,0)])
        self.encloses(r[0].val, -2.0)
        _, out, _ = OPS.accepted_update(
            matrix_identity(18), z, FULL._H_S(), matrix_identity(3), r)
        self.encloses(out[12].val, 4.0)

    def test_first_S_uses_estimator_mean_despite_large_truth_error(self):
        xhat = [Interval.point(0.0) for _ in range(18)]
        err = [Interval.point(0.0) for _ in range(18)]
        err[12] = Interval.point(300.0)
        _, out, _, cell, _ = FULL._measurement_branch_hull(
            matrix_identity(18), err, err[:3], FULL._H_S(), matrix_identity(3),
            [-xhat[12+i] for i in range(3)], allow_rejected=False)
        self.encloses(out[12], 300.0)
        self.assertLess(abs(cell["dx"][12].hi), 1e-300)

    def test_mean_plus_error_preserves_truth(self):
        xhat = [Interval.point(0.0) for _ in range(18)]
        xhat[12] = Interval.point(2.0)
        err = [Interval.point(0.0) for _ in range(18)]
        err[12] = Interval.point(3.0)
        _, ep, _, cell, _ = FULL._measurement_branch_hull(
            matrix_identity(18), err, err[:3], FULL._H_S(), matrix_identity(3),
            [-xhat[12+i] for i in range(3)], allow_rejected=False)
        xp = FULL._update_estimator_mean(xhat, cell, allow_rejected=False)
        self.encloses(xp[12], 1.0)
        self.encloses(ep[12], 4.0)
        self.encloses(xp[12]+ep[12], 5.0)

    def test_nominal_force_comes_from_nominal_aw(self):
        xhat = [Interval.point(0.0) for _ in range(18)]
        xhat[15] = Interval.point(20.0)
        f = FULL._predicted_force_upper(xhat, domain())
        self.assertGreaterEqual(f, 20.0+domain()["startup"]["gravity_mps2"])
        self.assertGreater(f, domain()["normal_live"]["specific_force_norm_upper_mps2"])
        self.assertGreaterEqual(FULL._H_acc(domain(), f)[0][1].hi, f)

    def test_prediction_consumes_current_gyro_bias_error(self):
        c = [Interval.point(0.0) for _ in range(3)]
        small = FULL._predict_c(c, matrix_identity(3), domain(), 0.005,
                                [Interval.point(0.01), *c[:2]])
        large = FULL._predict_c(c, matrix_identity(3), domain(), 0.005,
                                [Interval.point(0.2), *c[:2]])
        self.assertGreater(large[0].hi, 10*small[0].hi)

    def test_effective_acceleration_uses_vector_norm(self):
        e = [Interval.point(0.0) for _ in range(18)]
        for i in FULL.AW:
            e[i] = Interval.point(1.0)
        c = [Interval.point(0.0) for _ in range(3)]
        _, _, eta = FULL._acc_residual(e, c, domain(), 0.3, 0.0)
        required = FULL.VEFF.accel_latent_cross_gain_upper(0.3)*math.sqrt(3.0)
        self.assertGreaterEqual(eta, required)

    def test_horizontal_RS_squares_source_std_factor(self):
        factors = P3.source_rs_axis_std_factors()
        xy = struct.unpack("!f", struct.pack("!f", 0.72))[0]
        self.assertEqual(factors, [xy, xy, 1.0])
        rs = Interval.point(0.15)
        R = FULL._R_S({"R_S_filter_std": rs})
        self.encloses(R[0][0], (0.15*xy)**2)
        self.encloses(R[2][2], 0.15**2)
        self.assertLess(R[0][0].hi, R[2][2].lo)
        floor = P3.rs_variance_lower(rs, {"R_S_axis_std_factors": factors})
        self.assertLessEqual(floor, (0.15*xy)**2)
        self.assertGreater(floor, (0.15*xy)**2*(1-1e-12))

    def test_local_residual_jacobians_match_shipping_H(self):
        z = state()
        force = [Interval.point(x) for x in (1.0,2.0,-9.0)]
        mag = [Interval.point(x) for x in (20.0,5.0,40.0)]
        for residual, H in (
            (OPS.S_residual(z), FULL._H_S()),
            (OPS.accelerometer_residual(z, force), OPS.H_acc_canonical(force)),
            (OPS.magnetometer_residual(z, mag), OPS.H_mag_canonical(mag)),
        ):
            for i, row in enumerate(AD.jacobian(residual)):
                for j, a in enumerate(row):
                    self.encloses(a, H[i][j].lo)


CPP = r"""
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using F = Kalman3D_Wave_OU_III<double,true,false>;
using V = Eigen::Vector3d;
using Q = Eigen::Quaterniond;
F make_filter() {
    F f(V::Constant(0.02),V::Constant(0.001),V::Constant(0.25));
    f.qref=Q::Identity();
    f.xext.setZero();
    f.Pext.setIdentity();
    f.set_RS_noise(V::Ones());
    return f;
}
void print_c(const char* name, Q q) {
    std::cout << name;
    for(int i=0;i<3;++i) std::cout << " " << 2*q.vec()(i)/q.w();
    std::cout << "\n";
}
int main() {
    std::cout << std::setprecision(17);
    {
        SeaStateFusionFilter_OU_III<TrackerType::KALMANF> wrapper(false);
        // Call the actual source EMA/staging function. The time increment and
        // pending-consumption order are those in updateCore_.
        std::vector<int> staged, committed;
        for (int k=1; k<=64; ++k) {
            if (wrapper.online_tune_apply_pending_) committed.push_back(k);
            wrapper.online_tune_apply_pending_=false;
            wrapper.time_ += 0.005f;
            wrapper.adapt_mekf(0.005f,1.0f,0.1f,1.0f,2.0f);
            if (wrapper.online_tune_apply_pending_) staged.push_back(k);
        }
        std::cout << "staged_samples";
        for (int k:staged) std::cout << " " << k;
        std::cout << "\ncommitted_samples";
        for (int k:committed) std::cout << " " << k;
        std::cout << "\n";
    }
    {
        F f=make_filter();
        f.xext(12)=-1.0;
        f.Pext(0,0)=2.0;
        f.Pext(1,12)=f.Pext(12,1)=0.01;
        Q truth(1,0.2,0,0); truth.normalize();
        f.applyIntegralZeroPseudoMeas();
        // qref is WORLD->BODY. The declared error is R_true R_hat^T.
        print_c("S_c", truth*f.qref.conjugate());
        std::cout << "S_error " << -f.xext(12) << "\n";
        std::cout << "S_cov " << f.Pext(12,12) << "\n";
        std::cout << "S_cross_cov " << f.Pext(2,0) << "\n";
    }
    {
        F f=make_filter();
        f.xext(12)=2;
        f.applyIntegralZeroPseudoMeas();
        std::cout << "forced_S_error " << 5-f.xext(12) << "\n";
    }
    {
        F nominal=make_filter(), truth=make_filter();
        truth.qref=Q(1,0.2,-0.1,0.05).normalized();
        truth.xext.segment<3>(3)=V(0.007,-0.005,0.003);
        V omega(0.2,-0.3,0.1);
        nominal.time_update(omega,0.005);
        truth.time_update(omega,0.005);
        print_c("predict_c",truth.qref*nominal.qref.conjugate());
    }
}
"""


class CompiledFilterParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        inc = next((p for p in (ROOT/"third_party/eigen",Path("/usr/include/eigen3"))
                    if (p/"Eigen/Dense").is_file()), None)
        compiler = shutil.which("g++")
        if inc is None or compiler is None:
            msg = "C++ source parity requires g++ and Eigen"
            if os.environ.get("OU3_REQUIRE_CPP_PARITY") == "1":
                raise RuntimeError(msg)
            raise unittest.SkipTest(msg)
        with tempfile.TemporaryDirectory(prefix="ou3-source-parity-") as tmp:
            src=Path(tmp)/"golden.cpp"
            exe=Path(tmp)/"golden"
            src.write_text(CPP)
            subprocess.run([compiler,"-std=c++20","-O0",f"-I{ROOT/'src'}",
                            f"-I{inc}",str(src),"-o",str(exe)],check=True,timeout=180)
            output=subprocess.run([str(exe)],check=True,capture_output=True,
                                  text=True,timeout=30).stdout
        cls.golden={parts[0]:list(map(float,parts[1:]))
                    for line in output.splitlines() if (parts:=line.split())}
        print("CPP_SOURCE_PARITY="+json.dumps(cls.golden,sort_keys=True))

    def setUp(self):
        self.saved_n=FULL.N
        FULL.N=18

    def tearDown(self):
        FULL.N=self.saved_n

    def encloses(self, iv, value):
        # This tolerance only compares the C++ round-to-nearest diagnostic with
        # a real-arithmetic enclosure. It is never a theorem-promotion margin.
        tol=2e-14*max(1.0,abs(value))
        self.assertLessEqual(iv.lo,value+tol)
        self.assertGreaterEqual(iv.hi,value-tol)

    def test_AD_correction_matches_actual_filter(self):
        P,out,_=s_update()
        for i,value in enumerate(self.golden["S_c"]):
            self.encloses(out[i].val,value)
        self.encloses(out[12].val,self.golden["S_error"][0])
        self.encloses(P[12][12],self.golden["S_cov"][0])
        self.encloses(P[2][0],self.golden["S_cross_cov"][0])
        # Noncommuting rotation is essential: a collinear test misses the side.
        self.assertLess(self.golden["S_c"][2],-0.0009)

    def test_prefix_correction_matches_actual_filter(self):
        z=state({0:0.4,12:1.0})
        P=matrix_identity(18)
        P[0][0]=Interval.point(2.0)
        P[1][12]=P[12][1]=Interval.point(0.01)
        _,_,c,_,_=FULL._measurement_branch_hull(
            P,[x.val for x in z],[x.val for x in z[:3]],
            FULL._H_S(),matrix_identity(3),[x.val for x in OPS.S_residual(z)],
            allow_rejected=False)
        for i,value in enumerate(self.golden["S_c"]):
            self.encloses(c[i],value)

    def test_finite_bias_prediction_matches_actual_filter(self):
        d=domain()
        d["startup"]["effective_deterministic_gyro_transport_disturbance_upper_rad_s"]=0.0
        z=state({0:0.4,1:-0.2,2:0.1,3:0.007,4:-0.005,5:0.003})
        out=OPS.prediction(
            z,matrix_identity(18),matrix_identity(3),d,0.005,
            angular_rate_body=[Interval.point(x) for x in (0.2,-0.3,0.1)])
        for i,value in enumerate(self.golden["predict_c"]):
            self.encloses(out[i].val,value)

    def test_tuner_staging_clock_matches_actual_filter(self):
        now, last, pending = 0.0, 0.0, False
        staged, committed = [], []
        commit = PATH._constants()["commit"]
        for k in range(1, 65):
            out = PATH.source_clock_step(now, last, pending, 0.005, commit_s=commit)
            if out["commit_previous_candidate_before_sample"]:
                committed.append(k)
            pending = out["stage_updated_candidate_for_next_sample"]
            if pending:
                staged.append(k)
            now, last = out["time_s"], out["last_stage_s"]
        self.assertEqual(staged, self.golden["staged_samples"])
        self.assertEqual(committed, self.golden["committed_samples"])
        self.assertEqual(staged[:1], [21])
        self.assertEqual(committed[:1], [22])

    def test_forced_S_matches_actual_filter(self):
        z=state({12:3.0})
        _,out,_=OPS.accepted_update(
            matrix_identity(18),z,FULL._H_S(),matrix_identity(3),
            OPS.S_residual(z,truth_S=[Interval.point(x) for x in (5,0,0)]))
        self.encloses(out[12].val,self.golden["forced_S_error"][0])


if __name__ == "__main__":
    unittest.main()
