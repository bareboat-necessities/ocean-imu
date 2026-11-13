#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Reduced tangent geometry: (SO(3) semidirect R^(3N)) semidirect se(3).
// The action on the SIX bias coordinates is Ad_(R,v). This is the
// semi-direct-bias (SD) geometry, not the full 18-dimensional TG-INS with
// an extra virtual position-input bias. The CARMA dynamics are unchanged;
// no autonomous INS-error theorem is asserted for this stochastic hybrid.
#include "lie/TwoFrameGroup.h"
#include "lie/SE3Jacobians.h"

namespace ocean_imu::lie {
template<typename T, int NW = 6, int VCOL = 1>
struct TangentBiasGroup : TwoFrameGroup<T,NW,2> {
    using Base = TwoFrameGroup<T,NW,2>;
    using typename Base::Vector3;
    using typename Base::Matrix3;
    using typename Base::Tangent;
    using typename Base::AdjointMat;
    using Base::R; using Base::X; using Base::B;
    using Base::OFF_PHI; using Base::OFF_BG; using Base::OFF_BA;
    using Base::world_offset;
    using Base::Base;
    TangentBiasGroup() = default;
    explicit TangentBiasGroup(const Base& x) : Base(x) {}
    static TangentBiasGroup Identity() { return TangentBiasGroup(); }

    TangentBiasGroup compose(const TangentBiasGroup& y) const {
        TangentBiasGroup z(Base::compose(y));
        z.B.col(1) -= y.R.transpose() * Vector3(y.X.col(VCOL)).cross(Vector3(B.col(0)));
        return z;
    }
    TangentBiasGroup operator*(const TangentBiasGroup& y) const { return compose(y); }
    TangentBiasGroup inverse() const {
        TangentBiasGroup z(Base::inverse());
        z.B.col(1) -= Vector3(X.col(VCOL)).cross(Vector3(R * B.col(0)));
        return z;
    }
    static TangentBiasGroup Exp_G(const Tangent& x) {
        TangentBiasGroup z(Base::Exp_G(x));
        const Vector3 p = x.template segment<3>(OFF_PHI);
        const Vector3 v = x.template segment<3>(world_offset(VCOL));
        const Vector3 bg = x.template segment<3>(OFF_BG);
        const Vector3 jbg = left_jacobian<T>(p) * bg;
        const Matrix3 Ql = right_jacobian_cross<T>(Vector3(-p),Vector3(-v));
        z.B.col(1) += z.R.transpose() * (Ql * bg - Vector3(z.X.col(VCOL)).cross(jbg));
        return z;
    }
    static Tangent Log_G(const TangentBiasGroup& z) {
        Tangent x = Base::Log_G(z);
        const Vector3 p = x.template segment<3>(OFF_PHI);
        const Vector3 v = x.template segment<3>(world_offset(VCOL));
        const Vector3 bg = x.template segment<3>(OFF_BG);
        const Matrix3 Ql = right_jacobian_cross<T>(Vector3(-p),Vector3(-v));
        x.template segment<3>(OFF_BA) = inv_left_jacobian<T>(p) *
            (z.R * z.B.col(1) + Vector3(z.X.col(VCOL)).cross(Vector3(z.R*z.B.col(0))) - Ql*bg);
        return x;
    }
    TangentBiasGroup retract(const Tangent& x) const { return Exp_G(x).compose(*this); }
    Tangent local(const TangentBiasGroup& ref) const { return Log_G(compose(ref.inverse())); }
    AdjointMat Adjoint() const {
        AdjointMat A=Base::Adjoint();
        const Vector3 g=R*B.col(0);
        const Vector3 a=R*B.col(1)+Vector3(X.col(VCOL)).cross(g);
        const Matrix3 Vs=skew<T>(Vector3(X.col(VCOL))),Gs=skew<T>(g);
        A.template block<3,3>(OFF_BA,OFF_PHI)=(skew<T>(a)+Gs*Vs)*R;
        A.template block<3,3>(OFF_BA,OFF_BG)=Vs*R;
        A.template block<3,3>(OFF_BA,world_offset(VCOL))=Gs*R;
        return A;
    }

    // Differential coordinate shear from the two-frame error to SD error.
    static Tangent from_two_frame(const Tangent& x, const Vector3& v) {
        Tangent y=x;
        y.template segment<3>(OFF_BA) += v.cross(Vector3(x.template segment<3>(OFF_BG)));
        return y;
    }
    static AdjointMat ad(const Tangent& x) {
        AdjointMat A=AdjointMat::Zero();
        const Matrix3 P=skew<T>(Vector3(x.template head<3>()));
        for(int o=0;o<Base::NX;o+=3) A.template block<3,3>(o,o)=P;
        for(int o=3;o<Base::NX;o+=3)
            A.template block<3,3>(o,0)=skew<T>(Vector3(x.template segment<3>(o)));
        A.template block<3,3>(OFF_BA,world_offset(VCOL))=skew<T>(Vector3(x.template segment<3>(OFF_BG)));
        A.template block<3,3>(OFF_BA,OFF_BG)=skew<T>(Vector3(x.template segment<3>(world_offset(VCOL))));
        return A;
    }
    // The extra two-path term of J_r in row (bias_a), column (attitude).
    // Only 3x3 products, no dense matrix exponential or covariance temporary.
    static Matrix3 reset_extra(const Tangent& x) {
        const Matrix3 D=-skew<T>(Vector3(x.template head<3>()));
        const Matrix3 V=-skew<T>(Vector3(x.template segment<3>(world_offset(VCOL))));
        const Matrix3 G=-skew<T>(Vector3(x.template segment<3>(OFF_BG)));
        Matrix3 dp=Matrix3::Identity(),vp=Matrix3::Zero(),gp=Matrix3::Zero();
        Matrix3 zp=Matrix3::Zero(),sum=Matrix3::Zero();
        T coefficient=T(1);
        // At |phi| <= pi the remainder is below double roundoff. Translation
        // enters at most quadratically, so large wave units do not change
        // the rotational series' convergence radius.
        for(int n=1;n<=32;++n) {
            zp=(D*zp+G*vp+V*gp).eval();
            vp=(D*vp+V*dp).eval();
            gp=(D*gp+G*dp).eval();
            dp=(D*dp).eval();
            coefficient/=T(n+1);
            sum+=coefficient*zp;
        }
        return sum;
    }
    // In-place congruence by I + E_ba [C_phi C_bg C_v].
    template<class Mat>
    static void bias_row_congruence(Mat& P,const Matrix3& Cp,const Matrix3& Cg,const Matrix3& Cv) {
        constexpr int ov=world_offset(VCOL);
        for(int c=0;c<Base::NX;++c) {
            const Vector3 d=Cp*P.template block<3,1>(0,c)+Cg*P.template block<3,1>(OFF_BG,c)+Cv*P.template block<3,1>(ov,c);
            P.template block<3,1>(OFF_BA,c)+=d;
        }
        for(int r=0;r<Base::NX;++r) {
            const Eigen::Matrix<T,1,3> d=P.template block<1,3>(r,0)*Cp.transpose()+P.template block<1,3>(r,OFF_BG)*Cg.transpose()+P.template block<1,3>(r,ov)*Cv.transpose();
            P.template block<1,3>(r,OFF_BA)+=d;
        }
    }
};
} // namespace ocean_imu::lie
