// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include "ahrs/FrameConversions.h"
#include "ahrs/Mahony_AHRS.h"
#include <cstdio>
#include <algorithm>
int main() {
    float worst=0, roundtrip=0;
    for(float heading:{0.0f,90.0f,180.0f,-90.0f}) {
        Mahony_AHRS<float> filter;filter.init(0.2f,0.0f);
        auto truth=[&](float t){return quat_from_euler(8*std::sin(.61f*t),6*std::sin(.83f*t),heading+3*std::sin(.31415927f*t));};
        auto q=truth(0);const Eigen::Matrix3f D=Eigen::Vector3f(1,-1,-1).asDiagonal();
        Eigen::Quaternionf start(D*q.toRotationMatrix()*D);
        filter.q0=start.w();filter.q1=start.x();filter.q2=start.y();filter.q3=start.z();
        for(int k=1;k<=24000;++k) {
            const auto next=truth(k*.005f);Eigen::AngleAxisf delta(q.conjugate()*next);q=next;
            const Eigen::Quaternionf exact_nwu(D*q.toRotationMatrix()*D);
            const auto exact=mahony_nwu_to_bw_ned(exact_nwu);
            roundtrip=std::max(roundtrip,Eigen::AngleAxisf(exact.conjugate()*q).angle());
            const auto w=ned_to_mahony_nwu(delta.axis()*(delta.angle()/.005f));
            const auto a=ned_to_mahony_nwu(-9.80665f*(q.conjugate()*Eigen::Vector3f::UnitZ()));
            const auto m=ned_to_mahony_nwu(q.conjugate()*Eigen::Vector3f(20,0,47));
            float r,p,y;filter.updateMag(w.x(),w.y(),w.z(),a.x(),a.y(),a.z(),m.x(),m.y(),m.z(),&p,&r,&y,.005f);
            const auto got=mahony_nwu_to_bw_ned(Eigen::Quaternionf(filter.q0,filter.q1,filter.q2,filter.q3));
            const float err=Eigen::AngleAxisf(got.conjugate()*q).angle()*180/3.141592653589793f;
            worst=std::max(worst,err);
        }
    }
    std::printf("Mahony full-frame rotating-heading error %.6f deg\n",worst);
    // The exact basis identity is separate from Mahony's approximate inverse
    // square root and discrete feedback. A reversed 3-degree yaw oscillation
    // fails the tracking check by more than an order of magnitude.
    return roundtrip<1e-5f && worst<0.5f?0:1;
}
