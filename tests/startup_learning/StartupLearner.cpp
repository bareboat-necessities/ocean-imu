#include "StartupLearner.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace startup_learning {
Q rotation(const V& angle) {
    const double a = angle.norm();
    return a > 1e-12 ? Q(Eigen::AngleAxisd(a, angle / a))
                    : Q(1, angle.x()/2, angle.y()/2, angle.z()/2).normalized();
}
V rotation_log(const Q& q) {
    const Eigen::AngleAxisd a(q.normalized());
    return a.angle() * a.axis();
}
namespace {
struct Transport {
    std::vector<Q> q;
    Eigen::VectorXd residual;
    V field = V::Zero();
    std::size_t count = 0;
};
Transport transport(const std::vector<Sample>& samples, const V& bias, const V& mag_bias=V::Zero()) {
    Transport out;
    out.q.reserve(samples.size());
    Q q = Q::Identity();
    for (std::size_t i=0; i<samples.size(); ++i) {
        if (i) {
            const double dt = samples[i].time - samples[i-1].time;
            const V rate = 0.5 * (samples[i].gyro + samples[i-1].gyro) - bias;
            q = (q * rotation(dt * rate)).normalized();
        }
        out.q.push_back(q);
        if (samples[i].magnetic) { out.field += q * (samples[i].mag-mag_bias); ++out.count; }
    }
    if (out.count) out.field /= double(out.count);
    out.residual.resize(3 * out.count);
    std::size_t j = 0;
    for (std::size_t i=0; i<samples.size(); ++i)
        if (samples[i].magnetic) {
            out.residual.segment<3>(3*j) = out.q[i] * (samples[i].mag-mag_bias) - out.field;
            ++j;
        }
    return out;
}
V force_mean(const std::vector<Sample>& s, const Transport& tr, const V& ba, bool taper) {
    V sum = V::Zero(); double weight = 0;
    const double duration = s.back().time - s.front().time;
    for (std::size_t i=0; i<s.size(); ++i) {
        // Trapezoidal quadrature in physical time, also for irregular clocks.
        double dt = 0;
        if (i) dt += 0.5 * (s[i].time-s[i-1].time);
        if (i+1<s.size()) dt += 0.5 * (s[i+1].time-s[i].time);
        const double x = (s[i].time - s.front().time) / duration;
        const double w = dt * (taper ? 30*x*x*(1-x)*(1-x) : 1);
        sum += w * (tr.q[i] * (s[i].accel-ba)); weight += w;
    }
    return sum / weight;
}
V force_slope(const std::vector<Sample>& s, const Transport& tr, const V& ba) {
    V sum=V::Zero(); double norm=0;
    const double T=s.back().time-s.front().time;
    for(std::size_t i=0;i<s.size();++i) {
        double dt=0;
        if(i)dt+=0.5*(s[i].time-s[i-1].time);
        if(i+1<s.size())dt+=0.5*(s[i+1].time-s[i].time);
        const double x=(s[i].time-s.front().time)/T;
        const double w=30*x*x*(1-x)*(1-x)*dt;
        sum+=w*(x-.5)*T*(tr.q[i]*(s[i].accel-ba));
        norm+=w*(x-.5)*(x-.5)*T*T;
    }
    return sum/norm;
}
Eigen::VectorXd residual(const std::vector<Sample>& s,const Transport& tr,
                         const Config& c,bool joint) {
    Eigen::VectorXd r(tr.residual.size()+(joint?3:0));
    r.head(tr.residual.size())=tr.residual/c.mag_sigma;
    if(joint) {
        const double T=s.back().time-s.front().time;
        // A displacement nuisance enters this smooth slope as O(T^-3).
        // Preserve its finite-window uncertainty rather than treating 7000
        // accelerometer observations as independent gravity measurements.
        const double sigma=c.force_slope_sigma*std::pow(35/T,3);
        r.tail<3>()=force_slope(s,tr,c.accel_bias_prior)/sigma;
    }
    return r;
}
}
Learner::Learner(Config config) : config_(std::move(config)) {
    for(double value:{config_.gyro_bias_sigma,config_.accel_bias_sigma,config_.mag_sigma,
                      config_.mag_bias_sigma,config_.force_slope_sigma,config_.min_seconds,config_.max_seconds})
        if(!std::isfinite(value))throw std::invalid_argument("nonfinite startup configuration");
    if (!(config_.gyro_bias_sigma>0 && config_.accel_bias_sigma>0 &&
          config_.mag_sigma>0 && config_.mag_bias_sigma>0 && config_.force_slope_sigma>0 && config_.min_seconds>0 &&
          config_.max_seconds>=config_.min_seconds && config_.max_seconds<=35) ||
        !config_.gyro_bias_prior.allFinite() || !config_.accel_bias_prior.allFinite())
        throw std::invalid_argument("invalid startup configuration");
    samples_.reserve(7001);
}
void Learner::add(const Sample& s) {
    if (!std::isfinite(s.time) || !s.accel.allFinite() || !s.gyro.allFinite() ||
        (s.magnetic && (!s.mag.allFinite() || s.mag.norm()<1e-6)))
        throw std::invalid_argument("invalid startup sample");
    if (!samples_.empty()) {
        const double dt = s.time - samples_.back().time;
        if (!(dt>0 && dt<=0.1) || s.time-samples_.front().time>config_.max_seconds+1e-8)
            throw std::invalid_argument("startup clock gap, reversal or deadline exceeded");
    }
    samples_.push_back(s);
}
Handover Learner::estimate(const std::string& method) const {
    if (method!="uniform" && method!="tapered" && method!="magfit" && method!="joint" && method!="joint_cal")
        throw std::invalid_argument("unknown startup method");
    Handover out;
    if (samples_.size()<3 || samples_.back().time-samples_.front().time<1) return out;
    V bg = config_.gyro_bias_prior;
    V bm = V::Zero();
    const double precision = 1/(config_.gyro_bias_sigma*config_.gyro_bias_sigma);
    M bias_covariance = M::Identity()/precision;
    auto tr = transport(samples_,bg);
    if (tr.count<10) return out;
    if (method=="magfit" || method=="joint" || method=="joint_cal") {
        const bool joint=method!="magfit",calibrate=method=="joint_cal";
        const int size=calibrate?6:3;
        Eigen::VectorXd parameters=Eigen::VectorXd::Zero(size),prior=parameters;
        parameters.head<3>()=bg;prior.head<3>()=bg;
        Eigen::VectorXd precisions=Eigen::VectorXd::Constant(size,precision);
        if(calibrate)precisions.tail<3>().setConstant(1/(config_.mag_bias_sigma*config_.mag_bias_sigma));
        // Variable projection eliminates the unknown 3-D field. Repropagate
        // the entire gyro history for EACH bias candidate; never correct only
        // the terminal quaternion. The prior retains weak/unobserved axes.
        for (int iteration=0; iteration<5; ++iteration) {
            const Eigen::VectorXd r=residual(samples_,tr,config_,joint);
            Eigen::MatrixXd J(r.size(),size);
            constexpr double eps = 1e-5;
            for (int axis=0; axis<size; ++axis) {
                Eigen::VectorXd perturbed=parameters;perturbed[axis]+=eps;
                const V mb=calibrate?V(perturbed.tail<3>()):V::Zero();
                J.col(axis)=(residual(samples_,transport(samples_,perturbed.head<3>(),mb),config_,joint)-r)/eps;
            }
            Eigen::MatrixXd info=J.transpose()*J;
            info.diagonal()+=precisions;
            bias_covariance=info.ldlt().solve(Eigen::MatrixXd::Identity(size,size)).topLeftCorner<3,3>();
            const Eigen::VectorXd step = info.ldlt().solve(
                -J.transpose()*r - (precisions.array()*(parameters-prior).array()).matrix());
            const double old = r.squaredNorm() +
                               ((parameters-prior).array().square()*precisions.array()).sum();
            double scale=1;
            bool accepted=false;
            for (int search=0; search<8; ++search) {
                const Eigen::VectorXd proposal=parameters+scale*step;
                const V mb=calibrate?V(proposal.tail<3>()):V::Zero();
                auto trial = transport(samples_,proposal.head<3>(),mb);
                const double cost = residual(samples_,trial,config_,joint).squaredNorm() +
                    ((proposal-prior).array().square()*precisions.array()).sum();
                if (cost<=old) { parameters=proposal;tr=std::move(trial);accepted=true;break; }
                scale*=0.5;
            }
            if (!accepted || (scale*step).norm()<1e-8) break;
        }
        bg=parameters.head<3>();if(calibrate)bm=parameters.tail<3>();
    }
    const V f = force_mean(samples_,tr,config_.accel_bias_prior,method!="uniform");
    if (f.norm()<0.5*gravity) return out;
    const V down=-f.normalized();
    const V horizontal=tr.field-down*down.dot(tr.field);
    if (horizontal.norm()<1e-6) return out;
    const V north=horizontal.normalized(), east=down.cross(north);
    M frame; frame.row(0)=north.transpose(); frame.row(1)=east.transpose(); frame.row(2)=down.transpose();
    out.time=samples_.back().time;
    out.attitude=Q(frame*tr.q.back().toRotationMatrix()).normalized();
    out.field=frame*tr.field;
    out.down_body=out.attitude.conjugate()*V::UnitZ();
    out.gyro_bias=bg; out.accel_bias=config_.accel_bias_prior;
    out.mag_bias=bm;
    out.angular_rate=samples_.back().gyro-bg; out.last_raw_gyro=samples_.back().gyro;
    out.gyro_bias_cov=bias_covariance;
    out.magnetic_rms=std::sqrt(tr.residual.squaredNorm()/(3*tr.count));
    out.gravity_norm_error=std::abs(f.norm()-gravity);
    out.horizontal_field=horizontal.norm();
    const V uniform=force_mean(samples_,tr,config_.accel_bias_prior,false);
    out.gravity_sensitivity_rad=std::acos(std::clamp(down.dot(-uniform.normalized()),-1.,1.));
    out.valid=out.attitude.coeffs().allFinite() && out.field.allFinite();
    double last_mag_time=samples_.front().time;
    for(const auto& sample:samples_)if(sample.magnetic)last_mag_time=sample.time;
    // A low residual after inventing a bias far outside the supplied
    // calibration prior is model failure, not successful initialization.
    // This rejects magnetic steps that a constant-field fit explains by
    // spurious gyro drift / hard iron. It does not certify coherent offsets.
    const bool calibration_consistent=
        (bg-config_.gyro_bias_prior).norm()<=4*config_.gyro_bias_sigma &&
        bm.norm()<=4*config_.mag_bias_sigma;
    out.ready=out.valid && samples_.back().time-samples_.front().time>=config_.min_seconds &&
              out.horizontal_field>5 && out.magnetic_rms<2 && out.gravity_norm_error<0.5 &&
              out.gravity_sensitivity_rad<0.035 && samples_.back().time-last_mag_time<=0.2 &&
              calibration_consistent;
    return out;
}
} // namespace startup_learning
