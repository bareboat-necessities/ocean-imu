#define EIGEN_NON_ARDUINO
#define EIGEN_RUNTIME_NO_MALLOC
#include "imu_calibrate/MagTimingCal.h"
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include "util/MagneticInput.h"
#include "../common/MagTimingFixture.h"
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
using namespace ocean_imu::magnetic;
using namespace atoms3r_ical;
static void check(bool v,const char* why) { if(!v) { std::fprintf(stderr,"FAIL %s\n",why); std::exit(1); } }
static imu_cal::MagTimingFit replay(unsigned mode, unsigned lag=40000, uint32_t origin=1000000) {
    imu_cal::MagTimingCal cal;
    Q q=Q::Identity(); std::array<Q,64> poses{};
    V3 field(20,0,45), observed=field;
    uint32_t seq=0, stamp=origin;
    for(unsigned n=0;n<6000;++n) {
        const float t=n*.005f;
        V3 w(1.6f*std::sin(2.3f*t), 1.3f*std::sin(3.1f*t+.4f), 1.4f*std::cos(1.7f*t));
        if(mode==1) w=V3(0,0,1);
        if(mode==2) w.setZero();
        if(mode==5) w=V3(0,0,std::sin(3*t));
        if(n) q=(increment(w,.005f)*q).normalized();
        poses[n%64]=q;
        const uint32_t now=origin+n*5000;
        if(n>=24 && n%8==0) {
            unsigned delay=lag/5000;
            if(mode==3) delay+=unsigned((n/8)%5)-2; // +/-10 ms phase variation
            if(mode==4 && n>=3000) delay+=6; // validation disagrees by 30 ms
            observed=poses[(n-delay)%64]*field;
            if(mode==0 || mode==3) observed+=V3(.10f*std::sin(11*t),.13f*std::cos(17*t),.1f*std::sin(23*t));
            stamp=now; ++seq;
        }
        cal.add(now,w,seq!=0,seq,stamp,observed,n>=3000 ? 1 : 0);
    }
    const auto r=cal.fit();
    std::printf("mode=%u lag=%u ok=%d fit=%u width=%g allowance=%g rms=%g improve=%g pairs=%u/%u\n",
        mode,lag,r.ok,r.delay_us,r.fit_width_s,r.residual_time_s,r.validation_rms_rad,r.improvement,r.train_pairs,r.validation_pairs);
    return r;
}
struct KV {
    std::map<std::string,std::vector<uint8_t>> bytes;
    bool short_write=false;
    size_t getBytesLength(const char* k) { return bytes[k].size(); }
    size_t getBytes(const char* k,void* dst,size_t n) { if(bytes[k].size()!=n) return 0; std::memcpy(dst,bytes[k].data(),n); return n; }
    size_t putBytes(const char* k,const void* src,size_t n) { if(short_write) --n; const auto* p=static_cast<const uint8_t*>(src); bytes[k]={p,p+n}; return n; }
    bool remove(const char* k) { bytes.erase(k); return true; }
};
static void persistence() {
    ImuCalStoreT<KV> store;
    ImuCalBlobV5 b{};
    b.gyro_ok=1; b.gyro_coeff_crc=gyroCoeffCrc(b);
    b.mag_ok=1; b.mag_A[0]=b.mag_A[4]=b.mag_A[8]=1; b.mag_field_uT=50;
    b.sensor_id_lo=123; b.imu_type=3;
    check(gyroSetValid(b) && magSetValid(b),"spatial fixture");
    ImuCalBlobV5 loaded;
    check(!store.save(b),"v5 magnetic calibration requires timing");
    bindTimingFixture(b);
    check(store.saveVerified(b,loaded) && store.load(loaded),"complete v5 record round trip");
    check(loaded.mag_timing.delay_us==40000,"saved delay survives reload");
    const auto compatible=[&](uint32_t id,uint8_t type,uint8_t source,uint32_t config) {
        return magTimingCompatible(b,id,b.sensor_id_hi,type,source,config);
    };
    check(compatible(123,3,b.mag_source,b.mag_timing.configuration),"current timing can be preserved");
    check(!compatible(124,3,b.mag_source,b.mag_timing.configuration),"new device requires full calibration");
    check(!compatible(123,4,b.mag_source,b.mag_timing.configuration),"new IMU requires full calibration");
    check(!compatible(123,3,uint8_t(b.mag_source+1),b.mag_timing.configuration),"new source requires full calibration");
    check(!compatible(123,3,b.mag_source,b.mag_timing.configuration+1),"new config requires full calibration");
    auto changed=b; changed.mag_b[0]=1;
    check(!store.save(changed),"magnetic binding");
    changed=b; changed.gyro_b0[0]=.01f; changed.gyro_coeff_crc=gyroCoeffCrc(changed);
    check(!store.save(changed),"gyro binding");
    changed=b; changed.sensor_id_lo++;
    check(!store.save(changed),"device binding");
    changed=b; changed.accel_b0[0]=.1f;
    check(store.save(changed),"accel-only preserves timing");
    store.kv.bytes[store.kKeyV5][16]^=1;
    check(!store.load(loaded),"corruption rejected");
    changed=b; changed.mag_timing.residual_time_s=NAN;
    changed.mag_timing.crc=magTimingCrc(changed.mag_timing);
    check(!store.save(changed),"nonfinite timing rejected");
    store.kv.short_write=true;
    check(!store.save(b) && !store.load(loaded),"partial write rejected");
    store.kv.short_write=false;
    store.erase();
    // Old key is deliberately never read; no migration or fallback.
    const auto sealed=ImuCalStoreT<KV>::sealed_(b);
    store.kv.putBytes("blob_m5v4",&sealed,sizeof(sealed));
    check(!store.load(loaded),"old calibration key requires recalibration");
    changed=sealed; changed.version=4; changed.crc=computeBlobCrc(changed);
    store.kv.putBytes(store.kKeyV5,&changed,sizeof(changed));
    check(!store.load(loaded),"old version in new key rejected");

}
static void runtime() {
    Input input; input.reset(.00135f,.64f*M3::Identity());
    MagAcquisition timing;
    const V3 w(.5f,.3f,.7f), field(20,0,45), hard(8,-6,3);
    const uint32_t origin=0xffff0000u;
    Q expected=Q::Identity();
    for(unsigned n=0;n<=20;++n) {
        const uint32_t now=origin+n*5000;
        if(n==20) { timing.acquired(100,now,now,now+100); timing.delay_us=40000; timing.timing_allowance_s=.027f; }
        input.advance(now,w,1e-6f,timing,field+hard,hard);
        const Q proxy=increment(V3(.2f,0,0),-float(n)*.005f);
        input.capture(proxy,proxy,V3(float(n),0,-9.8f),w);
        if(n==12) expected=proxy;
    }
    V3 aligned; M3 covariance;
    check(input.currentField(hard,aligned,&covariance),"calibrated transport across clock wrap");
    check((aligned-increment(w,.04f)*field).norm()<1e-4f,"lag applied once before transport");
    check(input.observation().proxy_bw.angularDistance(expected)<1e-5f,"actual historical independent proxy");
    check(std::abs(input.observation().accel.x()-12)<1e-5f,"historical acceleration snapshot");
    check(input.observation().uncertainty.from_time_sigma>=.027f,"timing allowance propagated");
    check(covariance.trace()>1.92f,"uncertainty added to measurement covariance");
    timing.acquired(110,origin+110000,origin+110000,origin+110100); timing.delay_us=90000;
    input.advance(origin+110000,w,0,timing,field,hard);
    check(!input.currentField(hard,aligned),"invalid lag rejected");
}
int main() {
    Eigen::internal::set_is_malloc_allowed(false);
    for(unsigned lag:{20000u,40000u,60000u}) {
        const auto r=replay(0,lag); check(r.ok && std::abs(int(r.delay_us)-int(lag))<=2000,"known delay recovery");
    }
    check(!replay(1).ok,"constant-rate ambiguity");
    check(!replay(2).ok,"stillness ambiguity");
    check(!replay(5).ok,"single-axis excitation insufficient");
    check(!replay(4).ok,"held-out disagreement");
    const auto jitter=replay(3); check(jitter.ok && std::abs(int(jitter.delay_us)-40000)<=6000,"asynchronous jitter");
    check(replay(0,40000,0xffff0000u).ok,"calibration clock wrap");
    runtime();
    Eigen::internal::set_is_malloc_allowed(true);
    persistence();
    std::puts("Magnetic timing calibration PASS");
}
