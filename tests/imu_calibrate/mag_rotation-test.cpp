#define EIGEN_NON_ARDUINO
#define EIGEN_RUNTIME_NO_MALLOC
#include "util/MagneticRotation.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cstdio>
#include <cstdlib>
using namespace ocean_imu::magnetic;
static void check(bool ok, const char* text) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", text); std::exit(1); }
}
static void contracts() {
    RotationHistory<> h; RotationHistory<>::Rotation r;
    const V3 w(0,0,1), m(20,0,45);
    check(h.push(0,w)==Status::First, "zero timestamp is valid");
    for (uint32_t t=5000;t<=100000;t+=5000) h.push(t,w);
    check(h.between(60000,100000,r)==Status::Ready, "40 ms covered");
    check((r.q*m-increment(w,.04f)*m).norm()<1e-5f, "yaw rotation sign");
    check(h.push(100000,w)==Status::Duplicate, "duplicate gyro");
    check(h.push(99000,w)==Status::OutOfOrder, "out of order gyro");
    check(h.between(80000,100001,r)==Status::Future, "no extrapolation");
    check(h.between(100000,80000,r)==Status::OutOfOrder, "negative interval");
    check(h.between(uint32_t(-30000),100000,r)==Status::TooOld, "stale interval");
    check(h.push(130001,w)==Status::Gap, "missing gyro interval");
    check(h.between(90000,130001,r)==Status::NoHistory, "gap cannot be bridged");
    check(h.push(135001,V3(NAN,0,0))==Status::Invalid, "nonfinite gyro clears chain");
    h.reset(); uint32_t now=UINT32_MAX-20000; h.push(now,w);
    for(int i=0;i<10;++i) { now+=5000; h.push(now,w); }
    check(h.between(now-40000,now,r)==Status::Ready, "micros wrap");
    check((r.q*m-increment(w,.04f)*m).norm()<1e-5f, "wrap rotation");
    check(h.between(now-37500,now-2500,r)==Status::Ready, "partial boundary intervals");
    check((r.q*m-increment(w,.035f)*m).norm()<1e-5f, "boundary rate integration");
    h.reset(); h.push(2000000,w); h.push(2005000,w,.25f); h.push(2010000,w,.001f);
    check(h.between(2000000,2010000,r)==Status::Ready && r.bias_variance==.25f,
          "interval keeps the largest bias variance");
    atoms3r_ical::MagAcquisition a; a.acquired(7,100,110,120);
    check(a.sequence==1 && a.frame_us==100 && a.read_end_us==120 && a.fresh(207) && !a.fresh(208),
          "host metadata keeps the stale-source policy");
}
static void endurance() {
    RotationHistory<> h; RotationHistory<>::Rotation r;
    uint32_t t=0; const V3 w(.4f,-.7f,.9f); h.push(t,w);
    for(int i=0;i<2000000;++i) {t+=5000; h.push(t,w);}
    for(int i=0;i<1000;++i) check(h.between(t-120000,t,r)==Status::Ready, "long-duration history retained");
    check(std::abs(r.q.norm()-1)<2e-7f && r.visited<=32, "bounded lookup and norm after two wraps");
}
int main() {
    Eigen::internal::set_is_malloc_allowed(false);
    contracts(); endurance();
    std::puts("magnetic rotation history contracts PASS (Eigen allocation disabled)");
}
