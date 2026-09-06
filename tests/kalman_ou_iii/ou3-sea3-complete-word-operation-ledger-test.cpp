// Non-promoting selected-direction operation ledger for the admitted complete-SEA3 point.
//
// This translation unit reuses the checked-in complete-word feasibility harness so
// startup, source history, H18/A21 preparation, magnetic cadence, and the exact
// worst-direction coordinates are identical to that diagnostic.  It does NOT
// introduce a replay, finite harmonic source, independent F/Q schedule, tuner box,
// or independent R_S schedule.
//
// The ledger is a tangent diagnostic around the nominal shipping trajectory.  It
// reconstructs the operation-level covariance/error maps after executing the real
// wrapper call, and therefore must not be promoted to P4.  In particular, the
// difference between its tangent rho and the separately measured nonlinear rho is
// an explicit remaining nonlinear-transport obligation.

#define main ou3_complete_word_feasibility_embedded_main
#include "ou3-sea3-complete-word-feasibility-test.cpp"
#undef main

#include <array>
#include <cstdlib>
#include <string_view>

namespace {

constexpr int NX = 21;
constexpr int OFF_V = 6;
constexpr int OFF_S = 12;
constexpr int OFF_AW = 15;
constexpr int OFF_BA = 18;
using Vec21 = Eigen::Matrix<float, NX, 1>;
using Mat21 = Eigen::Matrix<float, NX, NX>;
using Mat21x3 = Eigen::Matrix<float, NX, 3>;
using Mat3x21 = Eigen::Matrix<float, 3, NX>;
using Vec12 = Eigen::Matrix<float, 12, 1>;

Mat21 reset_transport(const Eigen::Vector3f& dtheta) {
    Mat21 G = Mat21::Identity();
    if (dtheta.allFinite() && dtheta.squaredNorm() > 0.0f) {
        Eigen::Matrix3f S;
        S << 0.0f, -dtheta.z(), dtheta.y(),
             dtheta.z(), 0.0f, -dtheta.x(),
             -dtheta.y(), dtheta.x(), 0.0f;
        G.block<3,3>(0,0) = Eigen::Matrix3f::Identity() + 0.5f * S;
    }
    return G;
}

Mat21 joseph_from_pct(const Mat21& P,
                      const Mat21x3& K,
                      const Eigen::Matrix3f& S,
                      const Mat21x3& PCt) {
    Mat21 out = P - K * PCt.transpose() - PCt * K.transpose()
                  + K * S * K.transpose();
    return 0.5f * (out + out.transpose());
}

template<typename DerivedA, typename DerivedB>
double relative_matrix_difference(const Eigen::MatrixBase<DerivedA>& a,
                                  const Eigen::MatrixBase<DerivedB>& b) {
    return static_cast<double>((a.derived() - b.derived()).norm()) /
           std::max(1.0e-12, static_cast<double>(b.norm()));
}

double info_energy(const Mat21& Pfull, const Vec21& xfull, int dim) {
    if (dim == 18) {
        Eigen::Matrix<float,18,18> P = Pfull.topLeftCorner<18,18>();
        P = 0.5f * (P + P.transpose());
        Eigen::LDLT<Eigen::Matrix<float,18,18>> ldlt(P);
        if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0f)
            return std::numeric_limits<double>::quiet_NaN();
        const Eigen::Matrix<float,18,1> x = xfull.head<18>();
        return static_cast<double>(x.dot(ldlt.solve(x)));
    }
    Mat21 P = 0.5f * (Pfull + Pfull.transpose());
    Eigen::LDLT<Mat21> ldlt(P);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0f)
        return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(xfull.dot(ldlt.solve(xfull)));
}

Vec21 parse_direction(const std::string& raw, int dim) {
    Vec21 d = Vec21::Zero();
    std::stringstream ss(raw);
    std::string token;
    int i = 0;
    while (std::getline(ss, token, ',')) {
        if (i >= dim) throw std::runtime_error("too many direction entries");
        std::size_t used = 0;
        const float x = std::stof(token, &used);
        if (used != token.size() || !std::isfinite(x))
            throw std::runtime_error("invalid direction entry");
        d(i++) = x;
    }
    if (i != dim) throw std::runtime_error("direction dimension mismatch");
    return d;
}

struct EventAgg {
    int count = 0;
    double sum_delta = 0.0;
    double min_delta = std::numeric_limits<double>::infinity();
    double max_delta = -std::numeric_limits<double>::infinity();
    double max_ratio = 0.0;

    void add(double dv, double ratio) {
        ++count;
        sum_delta += dv;
        min_delta = std::min(min_delta, dv);
        max_delta = std::max(max_delta, dv);
        max_ratio = std::max(max_ratio, ratio);
    }
};

struct Ledger {
    int dim = 18;
    Vec21 error = Vec21::Zero();
    double V0 = std::numeric_limits<double>::quiet_NaN();
    double V1 = std::numeric_limits<double>::quiet_NaN();
    EventAgg prediction;
    EventAgg floor;
    EventAgg Szero;
    EventAgg accel;
    EventAgg vector;
    EventAgg lift;
    int imu_samples = 0;
    int vector_due = 0;
    int vector_accepted = 0;
    int vector_rejected = 0;
    double reconstruction_residual_max = 0.0;
    std::array<double,3> rs_std_min = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double,3> rs_std_max = {0.0, 0.0, 0.0};
    int rs_schedule_changes = 0;
    Eigen::Matrix3f last_RS = Eigen::Matrix3f::Constant(
        std::numeric_limits<float>::quiet_NaN());

    void observe_RS(const Eigen::Matrix3f& R) {
        if (!R.allFinite()) throw std::runtime_error("nonfinite committed R_S");
        for (int i = 0; i < 3; ++i) {
            if (!(R(i,i) > 0.0f)) throw std::runtime_error("nonpositive committed R_S diagonal");
            const double s = std::sqrt(static_cast<double>(R(i,i)));
            rs_std_min[static_cast<std::size_t>(i)] =
                std::min(rs_std_min[static_cast<std::size_t>(i)], s);
            rs_std_max[static_cast<std::size_t>(i)] =
                std::max(rs_std_max[static_cast<std::size_t>(i)], s);
        }
        if (last_RS.allFinite() && (R - last_RS).cwiseAbs().maxCoeff() > 0.0f)
            ++rs_schedule_changes;
        last_RS = R;
    }

    void event(EventAgg& a, const Mat21& Pbefore, const Mat21& Pafter,
               const Mat21& C) {
        const double vb = info_energy(Pbefore, error, dim);
        const Vec21 next = C * error;
        const double va = info_energy(Pafter, next, dim);
        if (!(std::isfinite(vb) && vb > 0.0 && std::isfinite(va) && va >= 0.0))
            throw std::runtime_error("invalid information energy in ledger event");
        a.add(va - vb, va / vb);
        error = next;
    }

    void dimension_lift(const Mat21& Pbefore, const Mat21& Pafter) {
        const double vb = info_energy(Pbefore, error, 18);
        const double va = info_energy(Pafter, error, 21);
        if (!(std::isfinite(vb) && vb > 0.0 && std::isfinite(va) && va >= 0.0))
            throw std::runtime_error("invalid H18-to-A21 lift energy");
        lift.add(va - vb, va / vb);
        dim = 21;
    }
};

void vector_event(Filter& f, Ledger& L, int word_k) {
    if (!due_mag(word_k)) return;
    ++L.vector_due;
    auto& m = f.mekf();
    const Mat21 Ppre = m.covariance_full();
    const bool active_pre = m.acc_bias_updates_enabled();
    f.updateMag(mag_truth());
    if (active_pre != m.acc_bias_updates_enabled())
        throw std::runtime_error("magnetic event changed fixed-dimensional mode");
    const auto& md = m.lastMagDiag();
    if (!md.accepted) {
        ++L.vector_rejected;
        return;
    }
    ++L.vector_accepted;
    const Mat21x3 PCt = m.PCt_scratch_;
    const Mat21x3 K = m.K_scratch_;
    Eigen::LDLT<Mat21> ldlt(Ppre);
    if (ldlt.info() != Eigen::Success)
        throw std::runtime_error("vector pre-covariance LDLT failed");
    const Mat21x3 Ht = ldlt.solve(PCt);
    const double lin_resid = relative_matrix_difference(Ppre * Ht, PCt);
    L.reconstruction_residual_max = std::max(L.reconstruction_residual_max, lin_resid);
    const Mat3x21 H = Ht.transpose();
    const Eigen::Vector3f dtheta = K.topRows<3>() * md.r;
    const Mat21 G = reset_transport(dtheta);
    const Mat21 C = G * (Mat21::Identity() - K * H);
    L.event(L.vector, Ppre, m.covariance_full(), C);
}

void imu_event(Filter& f, Ledger& L, int source_k, int word_k) {
    // #496 exposed the critical ordering bug: updateCore_ consumes a pending
    // online-tune commit before it builds the current Riccati/update packet.
    // Consume that exact shipping action once here, then snapshot F/Q/R_S.
    // The subsequent updateTime() sees no pending commit and is trajectory-equivalent.
    f.apply_pending_online_tune_();
    auto& m = f.mekf();

    const Mat21 P0 = m.covariance_full();
    Vec12 xlin0;
    xlin0.segment<3>(0) = m.get_velocity();
    xlin0.segment<3>(3) = m.get_position();
    xlin0.segment<3>(6) = m.get_integral_displacement();
    xlin0.segment<3>(9) = m.get_world_accel();

    const bool active_pre = m.acc_bias_updates_enabled();
    if ((L.dim == 18) != (!active_pre))
        throw std::runtime_error("ledger dimension disagrees with shipping bias mode");
    const bool floor_pending = m.aw_covariance_floor_pending_;
    const Eigen::Matrix3f floor_target = m.aw_covariance_floor_target_;
    const Eigen::Matrix3f R_S_pre = m.R_S;
    L.observe_RS(R_S_pre);
    const float pseudo_period_pre = m.pseudo_update_period_s_;
    const float tau_b_pre = m.tau_bacc_;
    const Eigen::Matrix3f Q_bacc_pre = m.Q_bacc_;
    float elapsed_copy = m.pseudo_update_elapsed_s_;
    const bool pseudo_due = ocean_imu::kalman::ou_detail::periodic_update_due(
        DT, pseudo_period_pre, elapsed_copy);

    f.updateTime(DT, gyro_zero(), acc_at(source_k), 35.0f);
    ++L.imu_samples;
    if (active_pre != m.acc_bias_updates_enabled())
        throw std::runtime_error("hybrid mode event inside fixed-dimensional IMU packet");

    Mat21 A = Mat21::Identity();
    A.block<6,6>(0,0) = m.F_AA_scratch_;
    A.block<12,12>(OFF_V,OFF_V) = m.F_LL_scratch_;
    const float tau_b = std::max(1.0e-3f, tau_b_pre);
    const float phi_b = active_pre ? std::exp(-DT / tau_b) : 1.0f;
    A.block<3,3>(OFF_BA,OFF_BA) = phi_b * Eigen::Matrix3f::Identity();

    Mat21 Q = Mat21::Zero();
    Q.block<6,6>(0,0) = m.Q_AA_scratch_;
    Q.block<12,12>(OFF_V,OFF_V) = m.Q_LL_scratch_;
    if (active_pre) {
        const float qd_scale = -0.5f * tau_b * std::expm1(-2.0f * DT / tau_b);
        Q.block<3,3>(OFF_BA,OFF_BA) = Q_bacc_pre * qd_scale;
    }

    Mat21 Pcur = A * P0 * A.transpose() + Q;
    Pcur = 0.5f * (Pcur + Pcur.transpose());
    L.event(L.prediction, P0, Pcur, A);

    if (floor_pending) {
        Eigen::Matrix3f Delta = floor_target - Pcur.block<3,3>(OFF_AW,OFF_AW);
        Delta = 0.5f * (Delta + Delta.transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(Delta);
        if (es.info() != Eigen::Success)
            throw std::runtime_error("a_w floor eigensolve failed");
        const Eigen::Vector3f evals = es.eigenvalues().cwiseMax(0.0f);
        const Eigen::Matrix3f DeltaPlus =
            es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
        Mat21 Pnext = Pcur;
        Pnext.block<3,3>(OFF_AW,OFF_AW) += DeltaPlus;
        L.event(L.floor, Pcur, Pnext, Mat21::Identity());
        Pcur = Pnext;
    }

    const Vec12 xlin_pred = m.F_LL_scratch_ * xlin0;
    if (pseudo_due) {
        Mat3x21 Hs = Mat3x21::Zero();
        Hs.block<3,3>(0,OFF_S) = Eigen::Matrix3f::Identity();
        Mat21x3 PCt = Pcur.block<NX,3>(0,OFF_S);
        if (!active_pre) PCt.block<3,3>(OFF_BA,0).setZero();
        const Eigen::Matrix3f S = Pcur.block<3,3>(OFF_S,OFF_S) + R_S_pre;
        Eigen::LDLT<Eigen::Matrix3f> ldlt(S);
        if (ldlt.info() != Eigen::Success)
            throw std::runtime_error("S=0 innovation LDLT failed");
        const Mat21x3 K = ldlt.solve(PCt.transpose()).transpose();
        const Eigen::Vector3f rS = -xlin_pred.segment<3>(6);
        const Eigen::Vector3f dtheta = K.topRows<3>() * rS;
        const Mat21 G = reset_transport(dtheta);
        const Mat21 C = G * (Mat21::Identity() - K * Hs);
        Mat21 Pnext = joseph_from_pct(Pcur, K, S, PCt);
        Pnext = G * Pnext * G.transpose();
        Pnext = 0.5f * (Pnext + Pnext.transpose());
        L.event(L.Szero, Pcur, Pnext, C);
        Pcur = Pnext;
    }

    const auto& ad = m.lastAccDiag();
    if (!ad.accepted)
        throw std::runtime_error("accelerometer update rejected inside admitted Normal-Live point word");
    const Mat21x3 PCt = m.PCt_scratch_;
    const Mat21x3 K = m.K_scratch_;
    Eigen::LDLT<Mat21> ldlt(Pcur);
    if (ldlt.info() != Eigen::Success)
        throw std::runtime_error("accelerometer pre-covariance LDLT failed");
    const Mat21x3 Ht = ldlt.solve(PCt);
    const double lin_resid = relative_matrix_difference(Pcur * Ht, PCt);
    const Mat3x21 H = Ht.transpose();
    const Eigen::Vector3f dtheta = K.topRows<3>() * ad.r;
    const Mat21 G = reset_transport(dtheta);
    const Mat21 C = G * (Mat21::Identity() - K * H);
    Mat21 Precon = joseph_from_pct(Pcur, K, m.S_scratch_, PCt);
    Precon = G * Precon * G.transpose();
    Precon = 0.5f * (Precon + Precon.transpose());
    const Mat21 Ppost = m.covariance_full();
    const double cov_resid = std::max(lin_resid, relative_matrix_difference(Precon, Ppost));
    L.reconstruction_residual_max = std::max(L.reconstruction_residual_max, cov_resid);
    if (!(std::isfinite(cov_resid) && cov_resid <= 2.0e-4))
        throw std::runtime_error("accelerometer covariance reconstruction mismatch");
    L.event(L.accel, Pcur, Ppost, C);

    vector_event(f, L, word_k);
}

void run_segment(Filter& f, Ledger& L, int source_start) {
    for (int k = 0; k < N; ++k)
        imu_event(f, L, source_start + k, k);
}

void emit_agg(std::ostream& os, const EventAgg& a) {
    os << '{'
       << "\"count\":" << a.count << ','
       << "\"sum_deltaV\":" << std::setprecision(17) << a.sum_delta << ','
       << "\"min_deltaV\":" << (a.count ? a.min_delta : 0.0) << ','
       << "\"max_deltaV\":" << (a.count ? a.max_delta : 0.0) << ','
       << "\"max_event_ratio\":" << a.max_ratio
       << '}';
}

void emit_ledger(std::ostream& os, const Ledger& L, const char* label,
                 double nonlinear_rho) {
    const double tangent_rho = L.V1 / L.V0;
    const double telescoped = L.prediction.sum_delta + L.floor.sum_delta +
        L.Szero.sum_delta + L.accel.sum_delta + L.vector.sum_delta + L.lift.sum_delta;
    os << '{'
       << "\"label\":\"" << label << "\",";
    os << "\"dimension_final\":" << L.dim << ',';
    os << "\"V0\":" << std::setprecision(17) << L.V0 << ',';
    os << "\"V1\":" << L.V1 << ',';
    os << "\"tangent_rho\":" << tangent_rho << ',';
    os << "\"nonlinear_rho_reference\":" << nonlinear_rho << ',';
    os << "\"tangent_minus_nonlinear_rho\":" << (tangent_rho - nonlinear_rho) << ',';
    os << "\"telescope_error\":" << (telescoped - (L.V1 - L.V0)) << ',';
    os << "\"imu_samples\":" << L.imu_samples << ',';
    os << "\"vector_due\":" << L.vector_due << ',';
    os << "\"vector_accepted\":" << L.vector_accepted << ',';
    os << "\"vector_rejected\":" << L.vector_rejected << ',';
    os << "\"rs_schedule_changes\":" << L.rs_schedule_changes << ',';
    os << "\"RS_std_min_xyz\":[" << L.rs_std_min[0] << ',' << L.rs_std_min[1] << ',' << L.rs_std_min[2] << "],";
    os << "\"RS_std_max_xyz\":[" << L.rs_std_max[0] << ',' << L.rs_std_max[1] << ',' << L.rs_std_max[2] << "],";
    os << "\"covariance_reconstruction_residual_max\":" << L.reconstruction_residual_max << ',';
    os << "\"prediction\":"; emit_agg(os, L.prediction); os << ',';
    os << "\"aw_floor\":"; emit_agg(os, L.floor); os << ',';
    os << "\"S_zero\":"; emit_agg(os, L.Szero); os << ',';
    os << "\"accelerometer\":"; emit_agg(os, L.accel); os << ',';
    os << "\"vector\":"; emit_agg(os, L.vector); os << ',';
    os << "\"H18_to_A21_lift\":"; emit_agg(os, L.lift);
    os << '}';
}

Ledger run_H(const Vec21& direction, bool continue_A) {
    auto f = prepare_mode('H');
    Ledger L;
    L.dim = 18;
    L.error = direction;
    L.V0 = info_energy(f->mekf().covariance_full(), L.error, 18);
    if (!(std::isfinite(L.V0) && L.V0 > 0.0))
        throw std::runtime_error("invalid H18 ledger initial energy");
    run_segment(*f, L, PRE);
    if (continue_A) {
        const Mat21 Pbefore = f->mekf().covariance_full();
        release_H_to_A(*f);
        const Mat21 Pafter = f->mekf().covariance_full();
        L.dimension_lift(Pbefore, Pafter);
        run_segment(*f, L, PRE + N);
    }
    L.V1 = info_energy(f->mekf().covariance_full(), L.error, L.dim);
    return L;
}

Ledger run_A(const Vec21& direction) {
    auto f = prepare_mode('A');
    Ledger L;
    L.dim = 21;
    L.error = direction;
    L.V0 = info_energy(f->mekf().covariance_full(), L.error, 21);
    if (!(std::isfinite(L.V0) && L.V0 > 0.0))
        throw std::runtime_error("invalid A21 ledger initial energy");
    run_segment(*f, L, PRE + N);
    L.V1 = info_energy(f->mekf().covariance_full(), L.error, 21);
    return L;
}

double nonlinear_connected_HA(const Vec21& direction) {
    auto nominal = prepare_mode('H');
    auto pert = prepare_mode('H');
    require_same_schedule(*nominal, *pert);
    const Eigen::MatrixXf P0 = covariance_mode(*nominal, 18);
    Eigen::VectorXf d18(18);
    d18 = direction.head<18>();
    inject(*pert, d18);
    const double V0 = energy(physical_increment(*pert, *nominal, 18), P0);
    if (!(std::isfinite(V0) && V0 > 0.0))
        throw std::runtime_error("invalid connected nonlinear initial energy");

    EventCounts cn, cp;
    run_H(*nominal, cn);
    run_H(*pert, cp);
    require_same_schedule(*nominal, *pert);
    release_H_to_A(*nominal);
    release_H_to_A(*pert);
    require_same_schedule(*nominal, *pert);
    run_A(*nominal, cn);
    run_A(*pert, cp);
    require_same_schedule(*nominal, *pert);

    const Eigen::MatrixXf P1 = covariance_mode(*nominal, 21);
    const double V1 = energy(physical_increment(*pert, *nominal, 21), P1);
    if (!(std::isfinite(V1) && V1 >= 0.0))
        throw std::runtime_error("invalid connected nonlinear final energy");
    return V1 / V0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            std::cerr << "usage: ou3-sea3-complete-word-operation-ledger-test <trace.txt> <H-dir> <H-rho> <A-dir> <A-rho>\n";
            return 2;
        }
        sfz.clear();
        std::ifstream in(argv[1]);
        if (!in) throw std::runtime_error("cannot open validated complete-SEA3 trace");
        float x;
        while (in >> x) {
            if (!std::isfinite(x)) throw std::runtime_error("nonfinite source input");
            sfz.push_back(x);
        }
        if (sfz.size() != static_cast<std::size_t>(TOTAL))
            throw std::runtime_error("complete-word trace length mismatch");

        const Vec21 hdir = parse_direction(argv[2], 18);
        const double hrho = std::stod(argv[3]);
        const Vec21 adir = parse_direction(argv[4], 21);
        const double arho = std::stod(argv[5]);
        if (!(std::isfinite(hrho) && hrho >= 0.0 && std::isfinite(arho) && arho >= 0.0))
            throw std::runtime_error("invalid nonlinear rho reference");

        const double harho = nonlinear_connected_HA(hdir);
        const Ledger H = run_H(hdir, false);
        const Ledger HA = run_H(hdir, true);
        const Ledger A = run_A(adir);

        std::cout << '{';
        std::cout << "\"qualification\":\"OU3_COMPLETE_SEA3_NONPROMOTING_OPERATION_LEDGER_V1\",";
        std::cout << "\"canonical_source\":\"COMPLETE_SEA3_NORMAL_LIVE_WORD\",";
        std::cout << "\"same_validated_continuum_member\":true,";
        std::cout << "\"pending_tune_consumed_at_shipping_boundary\":true,";
        std::cout << "\"actual_committed_per_axis_RS_observed\":true,";
        std::cout << "\"trajectory_replay_used\":false,";
        std::cout << "\"finite_harmonic_source_used\":false,";
        std::cout << "\"independent_F_Q_schedule_used\":false,";
        std::cout << "\"independent_RS_schedule_used\":false,";
        std::cout << "\"tangent_ledger_only\":true,";
        std::cout << "\"P4_promoted\":false,";
        std::cout << "\"H18\":"; emit_ledger(std::cout, H, "H18_worst_nonlinear_direction", hrho); std::cout << ',';
        std::cout << "\"H18_release_A21_connected\":"; emit_ledger(std::cout, HA, "H18_direction_connected_through_release_and_A21", harho); std::cout << ',';
        std::cout << "\"A21\":"; emit_ledger(std::cout, A, "A21_worst_nonlinear_direction", arho);
        std::cout << "}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "complete SEA3 operation-ledger failure: " << e.what() << '\n';
        return 1;
    }
}
