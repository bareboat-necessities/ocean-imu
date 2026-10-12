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
/**
 * Return the fit from two synthetic 15-second captures and check retry isolation.
 * Mode selects multi-axis motion (0), constant rate (1), stillness (2), phase
 * jitter (3), conflicting validation (4), single-axis rocking (5), field-axis
 * rotation (6), or noisy rocking with jitter (7) and additional pauses (8).
 * Lag and clock origin are microseconds; sensor strides count 5 ms ticks.
 */
static imu_cal::MagTimingFit replay(unsigned mode, unsigned lag=40000, uint32_t origin=1000000, unsigned mag_stride=8, unsigned imu_stride=1) {
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
        if(mode==6) w=field.normalized()*std::sin(3*t);
        if(mode==7) w=V3(.9f*std::sin(3.8f*t),0,0); // +/-14 degree rocking
        if(mode==8) w=V3(.9f*std::sin(3.8f*t),0,0);
        if(mode==8 && std::fmod(t,15.f)<10.f) w.setZero(); // long pauses
        if(n) q=(increment(w,.005f)*q).normalized();
        poses[n%64]=q;
        const uint32_t now=origin+n*5000;
        if(n>=24 && n%mag_stride==0) {
            unsigned delay=lag/5000;
            if(mode==3 || mode==7 || mode==8) delay+=unsigned((n/8)%5)-2; // +/-10 ms phase variation
            if(mode==4 && n>=3000) delay+=6; // validation disagrees by 30 ms
            observed=poses[(n-delay)%64]*field;
            if(mode==0 || mode==3 || mode>=6) observed+=V3(.10f*std::sin(11*t),.13f*std::cos(17*t),.1f*std::sin(23*t));
            stamp=now; ++seq;
        }
        if(n%imu_stride==0) cal.add(now,w,seq!=0,seq,stamp,observed,n>=3000 ? 1 : 0);
    }
    const auto r=cal.fit();
    check(!r.ok || cal.segmentQuality(0).ok,
          "accepted fit cannot be blocked by training progress");
    check(cal.segmentQuality(1).ok==r.ok && cal.segmentQuality(1).reason==r.reason,
          "validation progress uses the disjoint fit, including capture disagreement");
    std::printf("mode=%u lag=%u ok=%d reason=%s fit=%u width=%g allowance=%g rms=%g improve=%g pairs=%u/%u\n",
        mode,lag,r.ok,imu_cal::magTimingReason(r.reason),r.delay_us,r.fit_width_s,r.residual_time_s,r.validation_rms_rad,r.improvement,r.train_pairs,r.validation_pairs);
    const auto train_pairs=cal.pairs(0);
    cal.resetSegment(1);
    check(cal.pairs(0)==train_pairs && cal.pairs(1)==0 && !cal.fit().ok,
          "validation retry preserves training, cannot save without new validation");
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
    check(magTimingAbsent(b.mag_timing),"default timing record is absent");
    check(store.saveVerified(b,loaded) && store.load(loaded) && magTimingAbsent(loaded.mag_timing),
          "timing record is optional");
    auto partial=b; partial.mag_timing.delay_us=40000;
    check(!store.save(partial),"present timing record must be qualified");
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
/** Run delay, ambiguity and persistence regressions; exit on any failure. */
int main() {
    Eigen::internal::set_is_malloc_allowed(false);
    for(unsigned lag:{20000u,40000u,60000u}) {
        const auto r=replay(0,lag); check(r.ok && std::abs(int(r.delay_us)-int(lag))<=2000,"known delay recovery");
    }
    check(!replay(1).ok,"constant-rate ambiguity");
    check(!replay(2).ok,"stillness ambiguity");
    check(replay(5).ok,"varying single-axis motion identifies scalar delay");
    check(!replay(4).ok,"held-out disagreement");
    check(!replay(6).ok,"variable field-axis spin remains unobservable");
    for(unsigned mode:{7u,8u}) {
        const auto human=replay(mode,40000,1000000,8,2);
        check(human.ok && std::abs(int(human.delay_us)-40000)<=6000,
              "100 Hz moderate rocking with noise, phase jitter and pauses");
    }
    check(replay(0,40000,1000000,36,2).ok,"180 ms pairs need 260 ms history for all candidates");
    const auto jitter=replay(3); check(jitter.ok && std::abs(int(jitter.delay_us)-40000)<=6000,"asynchronous jitter");
    check(replay(0,40000,0xffff0000u).ok,"calibration clock wrap");
    Eigen::internal::set_is_malloc_allowed(true);
    persistence();
    std::puts("Magnetic timing calibration PASS");
}
