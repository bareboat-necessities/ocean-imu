/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Phase A of the right-invariant two-frame Lie-group filter: the group
    mathematics, on its own, before any estimator code exists.

    Two independent references are used, because a finite-difference check
    alone shares its assumptions with the code under test:

      1. central finite differences, for the Jacobians and the adjoint;
      2. Lie++ (group::SEn3<double,4>, third_party/Lie-plusplus), a separately
         derived implementation by the authors of the equivariant filtering
         literature this filter follows.

    Lie++ covers the world block only -- its SemiDirectBias is hardcoded to
    SE_2(3) + R^6 and cannot represent four world vectors plus six bias states.
    The bias block is therefore pinned by structural identities that a wrong
    Jacobian cannot satisfy: the one-parameter subgroup property, the exp/log
    round trip, and the exact adjoint relation.

    The whole battery runs at double (tight tolerances, the real check) and at
    float (loose tolerances, proving the firmware instantiation compiles and
    keeps its small-angle branches).
*/

#define EIGEN_NON_ARDUINO

#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

#include <groups/SEn3.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>

namespace lie = ocean_imu::lie;

namespace {

int failures = 0;

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

template<typename T>
bool close(T a, T b, T tol) {
    return std::abs(a - b) <= tol;
}

template<typename T, int R, int C>
bool mat_close(const Eigen::Matrix<T,R,C>& a, const Eigen::Matrix<T,R,C>& b, T tol) {
    return a.allFinite() && b.allFinite() && ((a - b).cwiseAbs().maxCoeff() <= tol);
}

template<typename T, int R, int C>
void report(const Eigen::Matrix<T,R,C>& a, const Eigen::Matrix<T,R,C>& b, const std::string& what) {
    std::cerr << "  " << what << " max|diff| = "
              << static_cast<double>((a - b).cwiseAbs().maxCoeff()) << '\n';
}

// ---------------------------------------------------------------------------
// Random sampling. Fixed seed: these are deterministic regressions, not a
// fuzzing campaign.
// ---------------------------------------------------------------------------

class Sampler {
public:
    explicit Sampler(unsigned seed) : rng_(seed) {}

    double uniform(double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(rng_);
    }

    template<typename T>
    Eigen::Matrix<T,3,1> vector3(double scale) {
        return Eigen::Matrix<T,3,1>(static_cast<T>(uniform(-scale, scale)),
                                    static_cast<T>(uniform(-scale, scale)),
                                    static_cast<T>(uniform(-scale, scale)));
    }

    // Rotation vector with a controlled magnitude, so tests can target the
    // near-identity, moderate, and near-pi regimes on purpose.
    template<typename T>
    Eigen::Matrix<T,3,1> rotation_vector(double angle) {
        Eigen::Matrix<T,3,1> axis = vector3<T>(1.0);
        if (axis.norm() < T(1e-6)) axis = Eigen::Matrix<T,3,1>(T(1), T(0), T(0));
        axis.normalize();
        return axis * static_cast<T>(angle);
    }

    template<typename Group>
    Group group_element(double angle, double scale) {
        using T = typename Group::Scalar;
        Group g;
        g.R = lie::Exp<T>(rotation_vector<T>(angle));
        for (int i = 0; i < Group::kWorldCols; ++i) g.X.col(i) = vector3<T>(scale);
        for (int j = 0; j < Group::kBiasCols; ++j) g.B.col(j) = vector3<T>(scale * 0.1);
        return g;
    }

    template<typename Group>
    typename Group::Tangent tangent(double angle, double scale) {
        using T = typename Group::Scalar;
        typename Group::Tangent xi;
        xi.setZero();
        xi.template segment<3>(Group::OFF_PHI) = rotation_vector<T>(angle);
        for (int i = 0; i < Group::kWorldCols; ++i) {
            xi.template segment<3>(Group::world_offset(i)) = vector3<T>(scale);
        }
        for (int j = 0; j < Group::kBiasCols; ++j) {
            xi.template segment<3>(Group::bias_offset(j)) = vector3<T>(scale * 0.1);
        }
        return xi;
    }

private:
    std::mt19937 rng_;
};

// Angles spanning the near-identity Taylor branch, the ordinary regime, and
// close to pi where Log's quaternion recovery has to stay well behaved.
const double kAngles[] = {0.0, 1e-9, 1e-6, 1e-3, 0.00999, 0.01001, 0.05, 0.5, 1.5, 3.0, 3.14};

// ---------------------------------------------------------------------------
// 1. SO(3) primitives
// ---------------------------------------------------------------------------

template<typename T>
void test_so3(Sampler& s, T tol, const std::string& tag) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using Vector3 = Eigen::Matrix<T,3,1>;

    for (double angle : kAngles) {
        const Vector3 phi = s.rotation_vector<T>(angle);
        const Matrix3 R = lie::Exp<T>(phi);
        const std::string at = " at angle " + std::to_string(angle);

        // Exp lands in SO(3).
        check(mat_close<T,3,3>(Matrix3(R * R.transpose()), Matrix3::Identity(), tol),
              tag + ": Exp is not orthogonal" + at);
        check(close<T>(R.determinant(), T(1), tol),
              tag + ": Exp determinant != 1" + at);

        // Exp/Log round trip. Near pi the axis sign is ambiguous, so compare
        // through the rotation rather than through the vector.
        const Vector3 phi_back = lie::Log<T>(R);
        check(mat_close<T,3,3>(lie::Exp<T>(phi_back), R, tol),
              tag + ": Exp/Log round trip" + at);
        if (angle < 3.0) {
            check(mat_close<T,3,1>(phi_back, phi, tol * T(10)),
                  tag + ": Log did not recover the rotation vector" + at);
        }

        // J_l J_l^-1 = I.
        const Matrix3 Jl  = lie::left_jacobian<T>(phi);
        const Matrix3 Jli = lie::inv_left_jacobian<T>(phi);
        check(mat_close<T,3,3>(Matrix3(Jl * Jli), Matrix3::Identity(), tol * T(10)),
              tag + ": J_l J_l^-1 != I" + at);

        // J_r(phi) = J_l(-phi) = Exp(-phi) J_l(phi). This identity is what the
        // bias block of the group exponential rests on.
        check(mat_close<T,3,3>(lie::left_jacobian<T>(Vector3(-phi)),
                               Matrix3(lie::Exp<T>(Vector3(-phi)) * Jl), tol * T(10)),
              tag + ": J_l(-phi) != Exp(-phi) J_l(phi)" + at);
        check(mat_close<T,3,3>(lie::right_jacobian<T>(phi),
                               lie::left_jacobian<T>(Vector3(-phi)), tol),
              tag + ": right_jacobian != J_l(-phi)" + at);
    }
}

// J_l(phi) = int_0^1 Exp(u phi) du, by Simpson. Double only: the quadrature
// error would swamp a float comparison.
void test_left_jacobian_quadrature(Sampler& s) {
    using T = double;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    for (double angle : kAngles) {
        const Eigen::Matrix<T,3,1> phi = s.rotation_vector<T>(angle);
        const int n = 2000;
        Matrix3 acc = Matrix3::Zero();
        for (int k = 0; k <= n; ++k) {
            const T u = static_cast<T>(k) / static_cast<T>(n);
            const T w = (k == 0 || k == n) ? T(1) : ((k % 2 == 1) ? T(4) : T(2));
            acc += w * lie::Exp<T>(Eigen::Matrix<T,3,1>(u * phi));
        }
        acc /= (T(3) * static_cast<T>(n));

        const Matrix3 Jl = lie::left_jacobian<T>(phi);
        if (!check(mat_close<T,3,3>(Jl, acc, 1e-10),
                   "J_l != integral of Exp(u phi) du at angle " + std::to_string(angle))) {
            report(Jl, acc, "J_l");
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Group axioms
// ---------------------------------------------------------------------------

template<typename Group>
void test_axioms(Sampler& s, typename Group::Scalar tol, const std::string& tag) {
    using T = typename Group::Scalar;
    using W = Eigen::Matrix<T,3,Group::kWorldCols>;
    using Bm = Eigen::Matrix<T,3,Group::kBiasCols>;

    auto same = [tol](const Group& a, const Group& b) {
        return mat_close<T,3,3>(a.R, b.R, tol) &&
               mat_close<T,3,Group::kWorldCols>(a.X, b.X, tol) &&
               mat_close<T,3,Group::kBiasCols>(a.B, b.B, tol);
    };

    for (int trial = 0; trial < 20; ++trial) {
        const Group a = s.group_element<Group>(s.uniform(0.0, 3.0), 5.0);
        const Group b = s.group_element<Group>(s.uniform(0.0, 3.0), 5.0);
        const Group c = s.group_element<Group>(s.uniform(0.0, 3.0), 5.0);
        const Group id = Group::Identity();

        check(same(a.compose(id), a), tag + ": X o I != X");
        check(same(id.compose(a), a), tag + ": I o X != X");

        // Inverse on both sides. The bias slot is what catches a law written
        // as a plain direct product.
        const Group ainv = a.inverse();
        const Group prods[2] = {a.compose(ainv), ainv.compose(a)};
        for (const Group& prod : prods) {
            check(mat_close<T,3,3>(prod.R, Eigen::Matrix<T,3,3>::Identity(), tol) &&
                  mat_close<T,3,Group::kWorldCols>(prod.X, W::Zero(), tol) &&
                  mat_close<T,3,Group::kBiasCols>(prod.B, Bm::Zero(), tol),
                  tag + ": X o X^-1 != I");
        }

        check(same(a.compose(b).compose(c), a.compose(b.compose(c))),
              tag + ": composition is not associative");
        check(same(a.compose(b).inverse(), b.inverse().compose(a.inverse())),
              tag + ": (ab)^-1 != b^-1 a^-1");
    }
}

// ---------------------------------------------------------------------------
// 3. Exp/Log, retraction, one-parameter subgroups
// ---------------------------------------------------------------------------

template<typename Group>
void test_exp_log(Sampler& s, typename Group::Scalar tol, const std::string& tag) {
    using T = typename Group::Scalar;

    for (double angle : kAngles) {
        if (angle > 3.0) continue;  // Log's axis is ambiguous at pi
        const std::string at = " at angle " + std::to_string(angle);
        const typename Group::Tangent xi = s.tangent<Group>(angle, 2.0);

        const Group g = Group::Exp_G(xi);
        const typename Group::Tangent back = Group::Log_G(g);
        if (!check(mat_close<T,Group::NX,1>(back, xi, tol * T(10)),
                   tag + ": Exp_G/Log_G round trip" + at)) {
            report(back, xi, "xi");
        }

        // Also check group -> tangent -> group, which catches an inverse
        // Jacobian that is wrong in a way the other direction absorbs.
        const Group h = s.group_element<Group>(angle, 3.0);
        const Group h_back = Group::Exp_G(Group::Log_G(h));
        check(mat_close<T,3,3>(h_back.R, h.R, tol) &&
              mat_close<T,3,Group::kWorldCols>(h_back.X, h.X, tol * T(10)) &&
              mat_close<T,3,Group::kBiasCols>(h_back.B, h.B, tol * T(10)),
              tag + ": Log_G/Exp_G round trip" + at);
    }

    typename Group::Tangent zero;
    zero.setZero();
    const Group e = Group::Exp_G(zero);
    check(mat_close<T,3,3>(e.R, Eigen::Matrix<T,3,3>::Identity(), tol) &&
          mat_close<T,3,Group::kWorldCols>(e.X, Eigen::Matrix<T,3,Group::kWorldCols>::Zero(), tol) &&
          mat_close<T,3,Group::kBiasCols>(e.B, Eigen::Matrix<T,3,Group::kBiasCols>::Zero(), tol),
          tag + ": Exp_G(0) != I");
}

/*
    One-parameter subgroup property:

        Exp_G(s1 xi) o Exp_G(s2 xi) = Exp_G((s1 + s2) xi)

    This is the structural test that pins the bias block. Lie++ cannot check it
    -- there is no oracle for four world vectors plus six body-frame biases --
    but a wrong J_l on the bias slot breaks it immediately, because the
    property only holds when the exponential really does integrate the ODE
    b'(s) = beta - [phi]x b(s) that the group law induces.
*/
template<typename Group>
void test_one_parameter_subgroup(Sampler& s, typename Group::Scalar tol, const std::string& tag) {
    using T = typename Group::Scalar;

    for (int trial = 0; trial < 20; ++trial) {
        const typename Group::Tangent xi = s.tangent<Group>(s.uniform(0.0, 1.2), 2.0);
        const T s1 = static_cast<T>(s.uniform(-1.0, 1.0));
        const T s2 = static_cast<T>(s.uniform(-1.0, 1.0));

        const Group lhs = Group::Exp_G(typename Group::Tangent(s1 * xi))
                              .compose(Group::Exp_G(typename Group::Tangent(s2 * xi)));
        const Group rhs = Group::Exp_G(typename Group::Tangent((s1 + s2) * xi));

        const bool ok = mat_close<T,3,3>(lhs.R, rhs.R, tol) &&
                        mat_close<T,3,Group::kWorldCols>(lhs.X, rhs.X, tol * T(10)) &&
                        mat_close<T,3,Group::kBiasCols>(lhs.B, rhs.B, tol * T(10));
        if (!check(ok, tag + ": Exp_G(s1 xi) Exp_G(s2 xi) != Exp_G((s1+s2) xi)")) {
            report(lhs.B, rhs.B, "bias block");
        }
    }
}

/*
    The retraction the filter actually applies, and its inverse -- "a
    correction followed by the inverse correction returns to where it started",
    which is what the update step depends on and the first thing to break if
    the injection lands on the wrong side.
*/
template<typename Group>
void test_retract_local(Sampler& s, typename Group::Scalar tol, const std::string& tag) {
    using T = typename Group::Scalar;
    using W = Eigen::Matrix<T,3,Group::kWorldCols>;
    using Bm = Eigen::Matrix<T,3,Group::kBiasCols>;

    for (int trial = 0; trial < 20; ++trial) {
        const Group g = s.group_element<Group>(s.uniform(0.0, 2.0), 4.0);
        const Group x = s.group_element<Group>(s.uniform(0.0, 2.0), 4.0);

        const typename Group::Tangent xi = x.local(g);
        const Group x_back = g.retract(xi);
        check(mat_close<T,3,3>(x_back.R, x.R, tol) &&
              mat_close<T,3,Group::kWorldCols>(x_back.X, x.X, tol * T(20)) &&
              mat_close<T,3,Group::kBiasCols>(x_back.B, x.B, tol * T(20)),
              tag + ": retract(local(x)) != x");

        const typename Group::Tangent eta = s.tangent<Group>(s.uniform(0.0, 1.0), 1.0);
        const typename Group::Tangent eta_back = g.retract(eta).local(g);
        check(mat_close<T,Group::NX,1>(eta_back, eta, tol * T(20)),
              tag + ": retract then local did not recover the correction");

        // Left injection, componentwise. This is the difference from the MEKF
        // pattern the redesign replaces, so assert the form explicitly.
        const Group e = Group::Exp_G(eta);
        const Group r = g.retract(eta);
        check(mat_close<T,3,3>(r.R, Eigen::Matrix<T,3,3>(e.R * g.R), tol),
              tag + ": retraction is not left multiplication on R");
        check(mat_close<T,3,Group::kWorldCols>(r.X, W(e.R * g.X + e.X), tol * T(10)),
              tag + ": world vectors did not rotate with the attitude correction");
        check(mat_close<T,3,Group::kBiasCols>(r.B, Bm(g.B + g.R.transpose() * e.B), tol * T(10)),
              tag + ": bias correction was not transported into the body frame");
    }
}

// ---------------------------------------------------------------------------
// 4. Adjoint
// ---------------------------------------------------------------------------

/*
    Exact relation, not a linearization:

        X Exp(xi) X^-1 = Exp(Ad_X xi)

    Checked at finite xi, so an adjoint that is only correct to first order
    fails here.
*/
template<typename Group>
void test_adjoint_exact(Sampler& s, typename Group::Scalar tol, const std::string& tag) {
    using T = typename Group::Scalar;

    for (int trial = 0; trial < 20; ++trial) {
        const Group X = s.group_element<Group>(s.uniform(0.0, 2.0), 4.0);
        const typename Group::Tangent xi = s.tangent<Group>(s.uniform(0.0, 1.0), 1.5);

        const Group lhs = X.compose(Group::Exp_G(xi)).compose(X.inverse());
        const Group rhs = Group::Exp_G(typename Group::Tangent(X.Adjoint() * xi));

        const bool ok = mat_close<T,3,3>(lhs.R, rhs.R, tol * T(10)) &&
                        mat_close<T,3,Group::kWorldCols>(lhs.X, rhs.X, tol * T(100)) &&
                        mat_close<T,3,Group::kBiasCols>(lhs.B, rhs.B, tol * T(100));
        if (!check(ok, tag + ": X Exp(xi) X^-1 != Exp(Ad_X xi)")) {
            report(lhs.X, rhs.X, "world block");
            report(lhs.B, rhs.B, "bias block");
        }
    }
}

// Central finite differences of the same relation. Double only.
void test_adjoint_finite_difference(Sampler& s) {
    using T = double;
    using Group = lie::TwoFrameGroup<T,4,2>;

    const T h = 1e-6;
    for (int trial = 0; trial < 10; ++trial) {
        const Group X = s.group_element<Group>(s.uniform(0.0, 2.0), 4.0);
        const Eigen::Matrix<T,Group::NX,Group::NX> Ad = X.Adjoint();

        Eigen::Matrix<T,Group::NX,Group::NX> numeric;
        for (int col = 0; col < Group::NX; ++col) {
            Group::Tangent e = Group::Tangent::Zero();
            e(col) = h;
            const Group plus  = X.compose(Group::Exp_G(e)).compose(X.inverse());
            e(col) = -h;
            const Group minus = X.compose(Group::Exp_G(e)).compose(X.inverse());
            numeric.col(col) = (Group::Log_G(plus) - Group::Log_G(minus)) / (T(2) * h);
        }

        if (!check(mat_close<T,Group::NX,Group::NX>(Ad, numeric, 1e-7),
                   "Adjoint does not match central finite differences")) {
            report(Ad, numeric, "Ad");
        }
    }
}

// ---------------------------------------------------------------------------
// 5. Lie++ cross-check of the world block
// ---------------------------------------------------------------------------

/*
    Our world block is SE_4(3): the (R, X_w) part of the two-frame group with
    the bias slot held at zero. Lie++'s group::SEn3<double,4> is that group,
    independently derived. Lie++'s tangent layout is [w, t_0, t_1, ...]
    contiguous; ours interleaves beta_g at offset 3. The world part is
    contiguous in both, so the mapping is a pair of segment copies.
*/
void test_liepp_world_block(Sampler& s) {
    using T = double;
    using Group = lie::TwoFrameGroup<T,4,2>;
    using Ref = group::SEn3<T,4>;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    auto idx = [](int i) { return static_cast<std::size_t>(i); };

    auto to_ref = [&idx](const Group& g) {
        typename Ref::IsometriesType t;
        for (int i = 0; i < 4; ++i) t[idx(i)] = g.X.col(i);
        return Ref(Matrix3(g.R), t);
    };
    auto tangent_to_ref = [](const Group::Tangent& xi) {
        Eigen::Matrix<T,15,1> u;
        u.segment<3>(0)  = xi.segment<3>(Group::OFF_PHI);
        u.segment<12>(3) = xi.segment<12>(Group::OFF_W);
        return u;
    };

    for (int trial = 0; trial < 25; ++trial) {
        Group a = s.group_element<Group>(s.uniform(0.0, 2.5), 6.0);
        Group b = s.group_element<Group>(s.uniform(0.0, 2.5), 6.0);
        a.B.setZero();
        b.B.setZero();

        // Composition.
        const Group ab = a.compose(b);
        const Ref ab_ref = to_ref(a) * to_ref(b);
        check(mat_close<T,3,3>(ab.R, Matrix3(ab_ref.R()), 1e-12),
              "Lie++ disagrees on composed rotation");
        for (int i = 0; i < 4; ++i) {
            check(mat_close<T,3,1>(Vector3(ab.X.col(i)), Vector3(ab_ref.t()[idx(i)]), 1e-12),
                  "Lie++ disagrees on composed world column " + std::to_string(i));
        }

        // Inverse.
        const Group ainv = a.inverse();
        const Ref ainv_ref = to_ref(a).inv();
        check(mat_close<T,3,3>(ainv.R, Matrix3(ainv_ref.R()), 1e-12),
              "Lie++ disagrees on inverse rotation");
        for (int i = 0; i < 4; ++i) {
            check(mat_close<T,3,1>(Vector3(ainv.X.col(i)), Vector3(ainv_ref.t()[idx(i)]), 1e-12),
                  "Lie++ disagrees on inverse world column " + std::to_string(i));
        }

        // Exponential.
        Group::Tangent xi = s.tangent<Group>(s.uniform(0.0, 2.0), 3.0);
        for (int j = 0; j < Group::kBiasCols; ++j) xi.segment<3>(Group::bias_offset(j)).setZero();
        const Group ex = Group::Exp_G(xi);
        const Ref ex_ref = Ref::exp(tangent_to_ref(xi));
        check(mat_close<T,3,3>(ex.R, Matrix3(ex_ref.R()), 1e-12),
              "Lie++ disagrees on Exp rotation");
        for (int i = 0; i < 4; ++i) {
            check(mat_close<T,3,1>(Vector3(ex.X.col(i)), Vector3(ex_ref.t()[idx(i)]), 1e-12),
                  "Lie++ disagrees on Exp world column " + std::to_string(i));
        }

        // Logarithm.
        const Group::Tangent lg = Group::Log_G(a);
        const Eigen::Matrix<T,15,1> lg_ref = Ref::log(to_ref(a));
        check(mat_close<T,3,1>(Vector3(lg.segment<3>(Group::OFF_PHI)),
                               Vector3(lg_ref.segment<3>(0)), 1e-12),
              "Lie++ disagrees on Log rotation part");
        check(mat_close<T,12,1>(Eigen::Matrix<T,12,1>(lg.segment<12>(Group::OFF_W)),
                                Eigen::Matrix<T,12,1>(lg_ref.segment<12>(3)), 1e-11),
              "Lie++ disagrees on Log world part");

        // Adjoint, world sub-blocks.
        const Eigen::Matrix<T,Group::NX,Group::NX> Ad = a.Adjoint();
        const Eigen::Matrix<T,15,15> Ad_ref = to_ref(a).Adjoint();
        check(mat_close<T,3,3>(Matrix3(Ad.block<3,3>(Group::OFF_PHI, Group::OFF_PHI)),
                               Matrix3(Ad_ref.block<3,3>(0,0)), 1e-12),
              "Lie++ disagrees on Adjoint attitude block");
        check(mat_close<T,12,3>(Eigen::Matrix<T,12,3>(Ad.block<12,3>(Group::OFF_W, Group::OFF_PHI)),
                                Eigen::Matrix<T,12,3>(Ad_ref.block<12,3>(3,0)), 1e-12),
              "Lie++ disagrees on Adjoint world/attitude coupling");
        check(mat_close<T,12,12>(Eigen::Matrix<T,12,12>(Ad.block<12,12>(Group::OFF_W, Group::OFF_W)),
                                 Eigen::Matrix<T,12,12>(Ad_ref.block<12,12>(3,3)), 1e-12),
              "Lie++ disagrees on Adjoint world block");
    }
}

// ---------------------------------------------------------------------------
// 6. Small-angle branch accuracy
// ---------------------------------------------------------------------------

/*
    Every small-angle branch in SO3Jacobians.h switches at 1e-2.

    Comparing f(seam - eps) against f(seam + eps) would NOT test this: those
    are different arguments, so the difference is dominated by the genuine
    slope of f, not by any branch mismatch. What has to hold is that *both*
    branches are accurate at the seam -- if each is within tol there, the
    function is continuous across it to within 2 tol.

    So evaluate against a closed form carried in long double, which has the
    headroom to absorb the cancellation in (1-cos t)/t^2 and 1 - (t/2)cot(t/2)
    at t = 1e-2 and still leave ~14 digits.
*/

using LD = long double;

Eigen::Matrix<LD,3,3> reference_exp(const Eigen::Matrix<LD,3,1>& phi) {
    using M = Eigen::Matrix<LD,3,3>;
    const LD t2 = phi.squaredNorm();
    const LD t = std::sqrt(t2);
    if (t == LD(0)) return M::Identity();
    const M W = lie::skew<LD>(phi);
    return M::Identity() + (std::sin(t) / t) * W + ((LD(1) - std::cos(t)) / t2) * (W * W);
}

Eigen::Matrix<LD,3,3> reference_left_jacobian(const Eigen::Matrix<LD,3,1>& phi) {
    using M = Eigen::Matrix<LD,3,3>;
    const LD t2 = phi.squaredNorm();
    const LD t = std::sqrt(t2);
    if (t == LD(0)) return M::Identity();
    const M W = lie::skew<LD>(phi);
    return M::Identity()
         + ((LD(1) - std::cos(t)) / t2) * W
         + ((t - std::sin(t)) / (t2 * t)) * (W * W);
}

Eigen::Matrix<LD,3,3> reference_inv_left_jacobian(const Eigen::Matrix<LD,3,1>& phi) {
    using M = Eigen::Matrix<LD,3,3>;
    const LD t2 = phi.squaredNorm();
    const LD t = std::sqrt(t2);
    if (t == LD(0)) return M::Identity();
    const M W = lie::skew<LD>(phi);
    const LD half = LD(0.5) * t;
    const LD c = (LD(1) - half * std::cos(half) / std::sin(half)) / t2;
    return M::Identity() - LD(0.5) * W + c * (W * W);
}

template<typename T>
void test_branch_seam(T tol, const std::string& tag) {
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    const Eigen::Matrix<LD,3,1> axis_ld =
        Eigen::Matrix<LD,3,1>(LD(0.3), LD(-0.5), LD(0.81)).normalized();
    const LD seam = static_cast<LD>(lie::kSmallAngle<T>);

    // Straddle the seam tightly, and sweep the decades either side so a branch
    // that is wrong away from the boundary is caught too.
    const LD angles[] = {
        LD(0), LD(1e-9), LD(1e-7), LD(1e-5), LD(1e-3),
        seam * LD(0.999), seam * LD(0.9999), seam,
        seam * LD(1.0001), seam * LD(1.001),
        LD(0.1), LD(0.7), LD(2.0), LD(3.0)
    };

    auto down = [](const Eigen::Matrix<LD,3,3>& m) { return Matrix3(m.template cast<T>()); };

    for (LD angle : angles) {
        const Eigen::Matrix<LD,3,1> phi_ld = axis_ld * angle;
        const Vector3 phi = phi_ld.cast<T>();
        const std::string at = " at angle " + std::to_string(static_cast<double>(angle));

        if (!check(mat_close<T,3,3>(lie::Exp<T>(phi), down(reference_exp(phi_ld)), tol),
                   tag + ": Exp disagrees with the closed form" + at)) {
            report(lie::Exp<T>(phi), down(reference_exp(phi_ld)), "Exp");
        }
        if (!check(mat_close<T,3,3>(lie::left_jacobian<T>(phi),
                                    down(reference_left_jacobian(phi_ld)), tol),
                   tag + ": left_jacobian disagrees with the closed form" + at)) {
            report(lie::left_jacobian<T>(phi), down(reference_left_jacobian(phi_ld)), "J_l");
        }
        // The cotangent form is genuinely singular at t = 0; there the Taylor
        // branch is the only truth and the reference cannot be evaluated.
        if (angle > LD(1e-8)) {
            if (!check(mat_close<T,3,3>(lie::inv_left_jacobian<T>(phi),
                                        down(reference_inv_left_jacobian(phi_ld)), tol),
                       tag + ": inv_left_jacobian disagrees with the closed form" + at)) {
                report(lie::inv_left_jacobian<T>(phi),
                       down(reference_inv_left_jacobian(phi_ld)), "J_l^-1");
            }
        }
        check(mat_close<T,3,3>(
                  Matrix3(lie::left_jacobian<T>(phi) * lie::inv_left_jacobian<T>(phi)),
                  Matrix3::Identity(), tol * T(10)),
              tag + ": J_l J_l^-1 != I" + at);
    }
}

// ---------------------------------------------------------------------------

template<typename T>
void run_all(T tol, const std::string& tag, bool with_fd_and_oracle) {
    Sampler s(20260806u);

    using Group  = lie::TwoFrameGroup<T,4,2>;   // full state: [v p S a_w], [b_g b_a]
    using GroupA = lie::TwoFrameGroup<T,4,1>;   // with_accel_bias = false

    static_assert(Group::NX == 21, "the full two-frame state must be 21-dimensional");
    static_assert(Group::OFF_PHI == 0 && Group::OFF_BG == 3 && Group::OFF_W == 6,
                  "tangent layout must match Kalman3D_Wave_OU_III");
    static_assert(Group::world_offset(0) == 6 && Group::world_offset(1) == 9 &&
                  Group::world_offset(2) == 12 && Group::world_offset(3) == 15,
                  "world columns must sit at the OU-III v/p/S/a_w offsets");
    static_assert(Group::bias_offset(0) == 3 && Group::bias_offset(1) == 18,
                  "bias columns must sit at the OU-III b_g/b_a offsets");
    static_assert(GroupA::NX == 18, "dropping the accelerometer bias must give 18 states");

    test_so3<T>(s, tol, tag);
    test_branch_seam<T>(tol, tag);
    test_axioms<Group>(s, tol, tag);
    test_axioms<GroupA>(s, tol, tag + " (no accel bias)");
    test_exp_log<Group>(s, tol, tag);
    test_one_parameter_subgroup<Group>(s, tol, tag);
    test_retract_local<Group>(s, tol, tag);
    test_adjoint_exact<Group>(s, tol, tag);

    if (with_fd_and_oracle) {
        test_left_jacobian_quadrature(s);
        test_adjoint_finite_difference(s);
        test_liepp_world_block(s);
    }
}

} // namespace

int main() {
    // Double is the real check. Float only has to compile, stay finite, and
    // keep its branches -- the firmware tolerance is set accordingly.
    run_all<double>(1e-11, "double", true);
    run_all<float>(2e-4f, "float", false);

    if (failures != 0) {
        std::cerr << failures << " Lie group check(s) failed\n";
        return 1;
    }
    std::cout << "lie_group-test: all checks passed\n";
    return 0;
}
