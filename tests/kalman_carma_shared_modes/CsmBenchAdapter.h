#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    carma3_bench adapter of the standalone shared-modes QMEKF
    (src/kalman_carma_shared_modes: KalmanCarmaSharedModes<float> behind
    the CARMA3 wrapper, SeaStateFusion_CarmaSharedModes) for
    kalman_carma_shared_modes-bench.  The engine, truth,
    sensors and metrics are carma3_bench's (tests/kalman_carma), so every
    filter of a row -- the CARMA3 baseline A41 and the filters below -- sees
    the same sea and the same IMU stream.

    Names: CSM[:noconf|:off][/option...]
      CSM          DEPLOYABLE default: estimated coupling + coupling confidence,
                   Evidence excitation (CSM@A: the Legacy excitation)
      CSM:noconf   confidence disabled (c = 1: the estimated coupling, unchecked)
      CSM:off      coupling off (the CARMA3 baseline A41 exactly)
    Robustness options (deliberately wrong coupling priors; the wrong prior is
    built from the true wave acceleration -- Sensors::a_wave_mag, the bench's
    ORACLE channel -- and enters as CouplingMode::External, so the filter
    must detect it from its own data):
      ext-true   the true band cross-spectrum (an oracle prior, for reference)
      ext-flip   x-z and y-z terms negated (wave direction reversed)
      ext-rot90  horizontal terms turned by 90 deg
      ext-ph45   x-z and y-z terms phase-shifted by 45 deg
      ext-amp2   x-z and y-z terms doubled (wrong amplitude; capped PSD by the mapping)
      fz<s>      freeze the estimated coupling at scenario time s ...
      keel       ... in the keel frame (it turns with the hull)
    Environment: CSM_TRACE=<path prefix> writes <prefix><label>.csv at 1 Hz
    (confidence, evidence, state, D^2, SNR, NIS ratio, f_op, sigma_z,
    applied coherence / phase, estimated direction).
*/

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "../kalman_carma/carma3/FullAdapters.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"

namespace carma3_bench {

namespace csm = ocean_imu::carma_shared_modes;

struct CsmSpec {
    bool confidence = true;
    bool off = false;
    std::string ext;              // "", true, flip, rot90, ph45, amp2
    float freeze_at_s = -1.0f;
    bool keel = false;
};

inline bool parse_csm_spec(const std::string& name, CsmSpec& sp) {
    if (name.rfind("CSM", 0) != 0) return false;
    const std::string rest = name.substr(3);
    std::vector<std::string> tok;
    size_t p = 0;
    while (p <= rest.size()) {
        const size_t q = rest.find('/', p);
        tok.push_back(rest.substr(p, q == std::string::npos ? std::string::npos : q - p));
        if (q == std::string::npos) break;
        p = q + 1;
    }
    if (tok[0] == ":noconf") sp.confidence = false;
    else if (tok[0] == ":off") sp.off = true;
    else if (!tok[0].empty()) return false;
    for (size_t i = 1; i < tok.size(); ++i) {
        const std::string& o = tok[i];
        if (o.rfind("ext-", 0) == 0) {
            sp.ext = o.substr(4);
            if (sp.ext != "true" && sp.ext != "flip" && sp.ext != "rot90" && sp.ext != "ph45" && sp.ext != "amp2")
                return false;
        } else if (o.rfind("fz", 0) == 0) sp.freeze_at_s = static_cast<float>(std::strtod(o.c_str() + 2, nullptr));
        else if (o == "keel") sp.keel = true;
        else return false;
    }
    return !(sp.off && !sp.ext.empty());
}

class FullCsm final : public IEstimator3 {
public:
    using Fusion = csm::SeaStateFusion_CarmaSharedModes;
    using Mekf = Fusion::Mekf;
    using H3 = ocean_imu::kalman::carma::Hermitian3<float>;

    FullCsm(std::string label, const CsmSpec& sp) : label_(std::move(label)), sp_(sp) {
        Fusion::Config fc;
        fill_config(fc, MAG_DELAY_SEC);
        Fusion::SharedModes c;
        c.mode = sp_.off ? csm::CouplingMode::Off
               : sp_.ext.empty() ? csm::CouplingMode::Estimated : csm::CouplingMode::External;
        c.confidence.enabled = sp_.confidence;
        c.confidence.accel_noise_std = ACC_SIGMA;
        f_.begin(fc, c);
        f_.raw().enableTuner(true);
        f_.raw().enableClamp(true);
        typename ocean_imu::kalman::carma::CarmaCrossSpectrumEstimator<float>::Config ec;
        ec.horizon_s = 600.0f;       // oracle prior (as the research ORACLE runs)
        ec.band_frac = c.band_frac;
        oracle_.configure(ec);
        if (sp_.keel) {
            auto& m = f_.raw().mekf().motion();
            auto cc = m.coupling().config();
            cc.frame = ocean_imu::kalman::carma::CarmaCouplingConfig::Frame::Keel;
            m.coupling().set_config(cc);
            // The facade re-applies its configuration with the gain each
            // commit; keep the keel frame there too.
            keel_ = true;
        }
        if (const char* tr = std::getenv("CSM_TRACE")) {
            std::string path = std::string(tr) + sanitize_(label_) + ".csv";
            trace_ = std::fopen(path.c_str(), "w");
            if (trace_)
                std::fprintf(trace_, "t,live,confidence,evidence,state,d2,snr,nis_ratio,coh_monitor,f_op,sigma_z,"
                                     "coh_xz,ph_xz,coh_yz,ph_yz,dir_to_deg,dir_unc_deg,shrinks\n");
        }
    }
    ~FullCsm() override { if (trace_) std::fclose(trace_); }

    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        temp_ = s.tempC;
        t_ += s.dt;
        mag_.tick(s.dt);
        f_.update(s.dt, s.gyr_ned, s.acc_ned, s.tempC);
        if (s.mag_tick) {
            f_.updateMag(s.mag_ned);
            mag_.after_mag(f_, f_.raw().mekf().lastMagDiag().accepted);
        }
        auto& m = f_.raw().mekf().motion();
        if (keel_) {
            auto cc = m.coupling().config();
            if (cc.frame != ocean_imu::kalman::carma::CarmaCouplingConfig::Frame::Keel) {
                cc.frame = ocean_imu::kalman::carma::CarmaCouplingConfig::Frame::Keel;
                m.coupling().set_config(cc);
            }
        }
        if (!sp_.ext.empty()) {
            oracle_.update(s.a_wave_mag.cast<float>(), s.dt, m.operating_frequency_hz());
            if (oracle_.elapsed_s() >= 30.0f) f_.raw().mekf().setExternalCoupling(wrong_(oracle_.cross()), 0.0f);
        }
        if (sp_.freeze_at_s >= 0.0f && t_ >= sp_.freeze_at_s && !frozen_) {
            m.coupling().freeze_now();
            frozen_ = true;
        }
        if (trace_ && t_ >= next_trace_) {
            next_trace_ += 1.0f;
            const auto& mk = f_.raw().mekf();
            const auto st = mk.coupling();
            const auto wd = mk.waveDirection();
            const double f_op = mk.motion().operating_frequency_hz(), sig_z = mk.motion().accel_std().z();
            std::fprintf(trace_, "%.2f,%d,%.4f,%.4f,%d,%.4g,%.4g,%.4f,%.4f,%.5f,%.5f,%.4f,%.2f,%.4f,%.2f,%.2f,%.2f,%lu\n",
                         static_cast<double>(t_), f_.isLive() ? 1 : 0, static_cast<double>(st.confidence.confidence),
                         static_cast<double>(st.confidence.evidence), static_cast<int>(st.confidence.state),
                         static_cast<double>(st.confidence.mismatch_d2), static_cast<double>(st.confidence.snr),
                         static_cast<double>(st.confidence.nis_ratio), static_cast<double>(st.confidence.monitor_coherence),
                         f_op, sig_z,
                         static_cast<double>(st.applied_coherence_xz), static_cast<double>(st.applied_phase_xz_deg),
                         static_cast<double>(st.applied_coherence_yz), static_cast<double>(st.applied_phase_yz_deg),
                         wd.valid ? static_cast<double>(wd.to_deg) : NAN,
                         wd.valid ? static_cast<double>(wd.uncertainty_deg) : NAN, st.shrink_events);
        }
    }
    Eigen::Matrix<double,9,9> native9_() const {
        const auto& m = f_.raw().mekf();
        const Eigen::MatrixXd P = m.covariance_full().template cast<double>();
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(9, Mekf::NX);
        for (int i = 0; i < 3; ++i) S(i, i) = 1.0;
        S.bottomRows(6) = m.pv_readout().template cast<double>();
        return S * P * S.transpose();
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = f_.raw().mekf();
        e.live = f_.isLive();
        e.p = m.get_position().template cast<double>();
        e.v = m.get_velocity().template cast<double>();
        e.a = m.get_world_accel().template cast<double>();
        const Eigen::Matrix<double,9,9> M9 = native9_();
        e.Ppv = M9.bottomRightCorner<6,6>();
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        e.has_ext = true;
        Eigen::Matrix<double,9,9> J = Eigen::Matrix<double,9,9>::Zero();
        J.topLeftCorner<3,3>() = e.q_bw.toRotationMatrix();
        J.block<3,3>(3,3) = -Eigen::Matrix3d::Identity();
        J.block<3,3>(6,6) = -Eigen::Matrix3d::Identity();
        e.P9 = J * M9 * J.transpose();
        e.b_g = m.gyroscope_bias().template cast<double>();
        e.b_a = m.get_acc_bias_body_at_temperature(temp_).template cast<double>();
        const auto& d = m.lastAccDiag();
        if (d.accepted) {
            Eigen::LLT<Eigen::Matrix3d> llt(d.S.template cast<double>());
            if (llt.info() == Eigen::Success) e.innov_white = llt.matrixL().solve(d.r.template cast<double>());
        }
        return e;
    }
    double native_nees9(const Truth9& tr) const override {
        const auto& m = f_.raw().mekf();
        const Eigen::Matrix3d R = m.quaternion_boat().template cast<double>().toRotationMatrix();
        Eigen::Matrix<double,9,1> err;
        const Eigen::AngleAxisd aa(Eigen::Matrix3d(tr.R_bw.transpose() * R));
        err.head<3>() = aa.angle() * aa.axis();
        err.segment<3>(3) = tr.p - m.get_position().template cast<double>();
        err.segment<3>(6) = tr.v - m.get_velocity().template cast<double>();
        Eigen::LDLT<Eigen::Matrix<double,9,9>> l(native9_());
        if (l.info() != Eigen::Success || !l.isPositive()) return NAN;
        return err.dot(l.solve(err));
    }
    double cond_corr() const override { return corr_condition(f_.raw().mekf().covariance_full()); }
    long mag_forwarded() const override { return mag_.forwarded; }
    long mag_applied() const override { return mag_.applied; }
    double op_fz_hz() const override {
        const auto& fr = f_.raw();
        return static_cast<double>(fr.getVerticalKappa()) * fr.mekf().vertical_model().base_zero_crossing_hz();
    }
    size_t memory_bytes() const override {
        return sizeof(f_) + sizeof(f_.raw().mekf());
    }
    int states() const override { return Mekf::NX; }
    void set_label(const std::string& l) override { label_ = l; }
    // <name>@<mode>[/opt...]: the CARMA excitation mode (default Evidence).
    bool set_options(const std::string& o) override { return apply_carma_excitation_options(f_.raw(), o); }
    bool excitation_diag(double (&d)[8]) const override { return carma_excitation_diag(f_.raw(), d); }
    bool coupling_diag(double (&d)[6]) const override {
        const H3& S = f_.raw().mekf().motion().excitation();
        const int pr[3][2] = {{0, 2}, {1, 2}, {0, 1}};
        for (int k = 0; k < 3; ++k) {
            const int i = pr[k][0], j = pr[k][1];
            const double den = static_cast<double>(S.re(i,i)) * S.re(j,j);
            const double re = S.re(i,j), im = S.im(i,j);
            d[2 * k] = den > 0 ? (re * re + im * im) / den : NAN;
            d[2 * k + 1] = std::atan2(im, re) * 180.0 / 3.14159265358979323846;
        }
        return true;
    }

private:
    static std::string sanitize_(std::string s) {
        for (char& ch : s) if (ch == '/' || ch == ':') ch = '_';
        return s;
    }
    H3 wrong_(H3 S) const {
        for (int i = 0; i < 2; ++i) {
            std::complex<float> v(S.re(i,2), S.im(i,2));
            if (sp_.ext == "flip") v = -v;
            else if (sp_.ext == "ph45") v *= std::polar(1.0f, 0.78539816f);
            else if (sp_.ext == "amp2") v *= 2.0f;
            S.re(i,2) = v.real(); S.im(i,2) = v.imag();
            S.re(2,i) = v.real(); S.im(2,i) = -v.imag();
        }
        if (sp_.ext == "rot90") {
            Eigen::Matrix3f R;
            R << 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;
            S.re = R * S.re * R.transpose();
            S.im = R * S.im * R.transpose();
        }
        return S;
    }

    std::string label_;
    CsmSpec sp_;
    Fusion f_;
    ocean_imu::kalman::carma::CarmaCrossSpectrumEstimator<float> oracle_;
    MagCounter mag_;
    float temp_ = 35.0f, t_ = 0.0f, next_trace_ = 0.0f;
    bool frozen_ = false, keel_ = false;
    std::FILE* trace_ = nullptr;
};

inline std::unique_ptr<IEstimator3> make_csm(const std::string& n, const CsmSpec& sp) {
    return std::make_unique<FullCsm>(n, sp);
}

}  // namespace carma3_bench
