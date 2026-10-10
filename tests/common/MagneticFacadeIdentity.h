#pragma once
// The three shipping startup facades must preserve their complete legacy
// startup/reference/refinement/continuous-hard-iron path at identity transport.
template <class Fusion, class CoreOf, class ProxyOf, class AttitudeOf, class Refined>
void magneticFacadeIdentity(CoreOf coreOf, ProxyOf proxyOf, AttitudeOf attitudeOf, Refined refined) {
    auto legacy=std::make_unique<Fusion>(), timed=std::make_unique<Fusion>();
    typename Fusion::Config cfg;
    legacy->begin(cfg); timed->begin(cfg);
    ocean_imu::magnetic::Observation observation;
    observation.covariance=cfg.sigma_m.array().square().matrix().asDiagonal();
    for(int i=1;i<=26000;++i) {
        const float t=float(i)*.005f;
        const V3 gyro(.03f*std::sin(t),.02f*std::cos(t),0), acc(0,0,-g_std);
        legacy->update(.005f,gyro,acc); timed->update(.005f,gyro,acc);
        if(i%8==0) {
            observation.proxy_bw=proxyOf(*timed); observation.attitude_bw=attitudeOf(*timed);
            observation.accel=acc; observation.gyro=gyro;
            legacy->updateMag(kField); timed->updateMagTransported(kField,observation);
        }
    }
    check(refined(*legacy)&&refined(*timed),"both facade paths refine magnetic reference");
    check((coreOf(*legacy).covariance_full()-coreOf(*timed).covariance_full()).norm()<2e-5f,
          "identity transport preserves complete startup/refinement/continuous path");
    check(coreOf(*legacy).quaternion().angularDistance(coreOf(*timed).quaternion())<1e-5f,
          "identity facade attitude parity");
}
