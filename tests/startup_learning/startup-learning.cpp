#include "ShippingBaseline.h"
#include "util/W3dSimCommon.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace sl = startup_learning;
namespace {
constexpr double rad = std::numbers::pi/180;
struct Truth { double time; sl::Q attitude; sl::V accel, gyro; };
std::vector<Truth> read_truth(const std::string& filename, double until) {
    std::ifstream in(filename); if (!in) throw std::runtime_error("cannot open input");
    std::string line; std::getline(in,line);
    std::vector<Truth> records;
    while (std::getline(in,line)) {
        std::istringstream row(line); double a[23]; char comma;
        for (int j=0;j<23;++j) {
            if (!(row>>a[j]) || !std::isfinite(a[j])) throw std::runtime_error("invalid CSV value");
            if(j<22 && (!(row>>comma)||comma!=',')) throw std::runtime_error("invalid CSV separator");
        }
        if(a[0]>until+0.001) break;
        if(!records.empty() && std::abs(a[0]-records.back().time-.005)>1e-8)
            throw std::runtime_error("expected a contiguous 200 Hz record");
        sl::Q q(a[19],a[20],a[21],a[22]);
        if (std::abs(q.norm()-1)>1e-5) throw std::runtime_error("invalid truth quaternion");
        const sl::M S=basis_zu_to_ned().cast<double>();
        Truth t{a[0], sl::Q(S*q.normalized().toRotationMatrix().transpose()*S),
                S*sl::V(a[10],a[11],a[12]), S*sl::V(a[13],a[14],a[15])};
        records.push_back(t);
    }
    return records;
}
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out; std::istringstream in(s); std::string x;
    while(std::getline(in,x,',')) out.push_back(x);
    if(out.empty()) throw std::invalid_argument("empty list");
    return out;
}
double angle(const sl::V& a, const sl::V& b) {
    return std::acos(std::clamp(a.normalized().dot(b.normalized()),-1.,1.))/rad;
}
double yaw(const sl::Q& q) {const auto R=q.toRotationMatrix();return std::atan2(R(1,0),R(0,0))/rad;}
void emit(std::ostream& out,const std::string& input,const std::string& profile,unsigned seed,
          double offset,const std::string& method,const std::string& schedule,
          const sl::Handover& h,const Truth& truth,const sl::V& bg,const sl::V& ba,double ms) {
    const sl::Q magnetic_from_geo(Eigen::AngleAxisd(-MagSim_WMM::default_declination_deg*rad,sl::V::UnitZ()));
    const sl::Q qtruth=magnetic_from_geo*truth.attitude;
    const sl::V B=magnetic_from_geo*MagSim_WMM::mag_world_aero().cast<double>();
    const sl::V be=h.gyro_bias-bg, ae=h.accel_bias-ba, rate=h.angular_rate-truth.gyro;
    out<<input<<','<<profile<<','<<seed<<','<<offset<<','<<method<<','<<schedule<<','
       <<h.time-offset<<','<<h.valid<<','<<h.ready<<','
       <<diffDeg(float(yaw(h.attitude)),float(yaw(qtruth)))<<','
       <<angle(h.down_body,qtruth.conjugate()*sl::V::UnitZ())<<','
       <<sl::rotation_log(qtruth.conjugate()*h.attitude).norm()/rad<<','
       <<angle(h.field,B)<<','<<h.field.norm()-B.norm()<<','<<(h.field-B).norm()<<','
       <<be.x()/rad<<','<<be.y()/rad<<','<<be.z()/rad<<','<<be.norm()/rad<<','
       <<ae.x()<<','<<ae.y()<<','<<ae.z()<<','<<ae.norm()<<','<<rate.norm()/rad<<','
       <<h.magnetic_rms<<','<<h.gravity_norm_error<<','<<h.gravity_sensitivity_rad/rad<<','<<ms;
    // Reviewable handover values, not only errors against hidden truth.
    out<<','<<h.time<<','<<h.attitude.w()<<','<<h.attitude.x()<<','<<h.attitude.y()<<','<<h.attitude.z();
    for(const auto& v:{h.field,h.down_body,h.gyro_bias,h.accel_bias,h.mag_bias,h.angular_rate,h.last_raw_gyro})
        for(int i=0;i<3;++i)out<<','<<v[i];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)out<<','<<h.gyro_bias_cov(i,j);
    out<<'\n';
}
}
int main(int argc,char** argv) try {
    std::string input,output,seed_text="1",offset_text="0",profile_text="calibrated";
    for(int i=1;i<argc;++i) {
        const std::string key=argv[i];
        if(i+1>=argc) throw std::invalid_argument("option requires a value");
        const std::string value=argv[++i];
        if(key=="--input")input=value;else if(key=="--output")output=value;
        else if(key=="--seeds")seed_text=value;else if(key=="--offsets")offset_text=value;
        else if(key=="--profiles")profile_text=value;else throw std::invalid_argument("unknown option: "+key);
    }
    if(input.empty()||output.empty())throw std::invalid_argument("--input and --output required");
    std::vector<double> offsets;for(const auto& x:split(offset_text)) {
        const double v=std::stod(x);if(!std::isfinite(v)||v<0)throw std::invalid_argument("invalid offset");
        offsets.push_back(v);
    }
    const auto records=read_truth(input,*std::max_element(offsets.begin(),offsets.end())+35);
    std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);out.open(output);
    out<<std::setprecision(10);
    out<<"input,profile,seed,offset,method,schedule,elapsed_s,valid,ready,yaw_error_deg,tilt_error_deg,attitude_error_deg,field_direction_error_deg,field_magnitude_error_uT,field_error_uT,bgx_error_dps,bgy_error_dps,bgz_error_dps,bg_error_dps,bax_error_mps2,bay_error_mps2,baz_error_mps2,ba_error_mps2,rate_error_dps,mag_residual_uT,gravity_norm_error_mps2,gravity_sensitivity_deg,solve_ms,epoch_s,qw,qx,qy,qz";
    for(const std::string name:{"field_uT","down_body","bg_rps","ba_mps2","bm_uT","omega_rps","last_gyro_rps"})
        for(const char axis:{'x','y','z'})out<<','<<name<<'_'<<axis;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)out<<",pbg_"<<i<<j;
    out<<'\n';
    std::size_t trials=0;
    for(const auto& profile:split(profile_text)) for(const auto& seed_string:split(seed_text)) for(double offset:offsets) {
        if(profile!="ideal"&&profile!="calibrated"&&profile!="residual"&&profile!="interference")
            throw std::invalid_argument("unknown noise profile");
        const unsigned seed=std::stoul(seed_string);
        const bool noise=profile!="ideal", residual=profile=="residual";
        auto an=make_imu_noise_model(noise?1.51e-3f*g_std:0,noise?8e-3f*g_std:0,noise?.0005f:0,w3d_expand_seed(seed,11),w3d_expand_seed(seed,12));
        auto gn=make_imu_noise_model(noise?.00157f:0,noise?float(.05*rad):0,noise?.00001f:0,w3d_expand_seed(seed,13),w3d_expand_seed(seed,14));
        auto mn=make_mag_noise_model(noise?.8f:0,residual?2.f:0,residual?.01f:0,residual?.015f:0,residual?.010f:0,residual?1.f:0,w3d_expand_seed(seed,15),w3d_expand_seed(seed,16));
        sl::Learner learner;
        sl::ShippingBaseline shipping;
        std::size_t count=0;
        for(const auto& truth:records) {
            if(truth.time<offset-1e-8) continue;
            const double t=truth.time-offset;if(t>35+1e-8)break;
            sl::Sample sample;
            sample.time=truth.time;
            sample.accel=apply_imu_noise(truth.accel.cast<float>(),an,.005f).cast<double>();
            sample.gyro=apply_imu_noise(truth.gyro.cast<float>(),gn,.005f).cast<double>();
            // Deterministic 25 Hz acquisition clock; no repeated held samples.
            sample.magnetic=(count%8==0);
            if(sample.magnetic) {
                const sl::V m=truth.attitude.conjugate()*MagSim_WMM::mag_world_aero().cast<double>();
                sample.mag=apply_mag_noise(m.cast<float>(),mn,.04f).cast<double>();
                if(profile=="interference"&&t>=15)sample.mag+=sl::V(10,0,0);
            }
            learner.add(sample);
            const sl::V bg=(gn.turn_on_bias+gn.bias_rw).cast<double>();
            const sl::V ba=(an.turn_on_bias+an.bias_rw).cast<double>();
            if(shipping.step(sample)) {
                emit(out,std::filesystem::path(input).filename(),profile,seed,offset,"shipping","handoff",shipping.handover(),truth,bg,ba,0);
                for(const std::string method:{"uniform","tapered","magfit","joint","joint_cal"}) {
                    const auto start=std::chrono::steady_clock::now();const auto h=learner.estimate(method);
                    const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
                    emit(out,std::filesystem::path(input).filename(),profile,seed,offset,method,"matched",h,truth,bg,ba,ms);
                }
            }
            if(count>=2000&&count%1000==0)for(const std::string method:{"uniform","tapered","magfit","joint","joint_cal"}) {
                const auto start=std::chrono::steady_clock::now();const auto h=learner.estimate(method);
                const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
                emit(out,std::filesystem::path(input).filename(),profile,seed,offset,method,"fixed",h,truth,bg,ba,ms);
            }
            ++count;
        }
        if(count!=7001)throw std::runtime_error("incomplete 35-second startup window");
        ++trials;
    }
    out.close();std::cerr<<std::filesystem::path(input).filename()<<": "<<trials<<" paired windows\n";
    return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
