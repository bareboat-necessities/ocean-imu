// Synthetic workload benchmark. Does not read/write calibration records.
#include <Arduino.h>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>
#include <esp_heap_caps.h>
#include <esp_timer.h>
#include <new>
#include "imu_calibrate/GyroCalFit.h"
#include "imu_calibrate/SensorNoise.h"
#include "kalman_common/MagneticInnovationGate.h"

using V=Eigen::Vector3f;
using M=Eigen::Matrix3f;
using Fit=imu_cal::GyroRotationFit<32768,24>;
static uint32_t random_word=17;
static float uniformNoise() {
  random_word=1664525u*random_word+1013904223u;
  return (float(random_word>>8)/16777216.f-.5f)*3.464101615f;
}
static void benchmark(void*) {
  Serial.printf("[BENCH] CPU_MHz=%u free_internal=%u free_PSRAM=%u\n",getCpuFrequencyMhz(),
    unsigned(heap_caps_get_free_size(MALLOC_CAP_INTERNAL)),unsigned(heap_caps_get_free_size(MALLOC_CAP_SPIRAM)));
  void* memory=heap_caps_malloc(sizeof(Fit),MALLOC_CAP_SPIRAM|MALLOC_CAP_8BIT);
  if(!memory){Serial.println("[BENCH] enable OPI PSRAM; no buffer available");vTaskDelete(nullptr);return;}
  Fit* fit=new(memory) Fit{};
  M S;S<<1.035f,.018f,-.012f,-.009f,.975f,.022f,.014f,-.017f,1.02f;
  const V bias(.006f,-.009f,.004f);uint32_t us=10000;
  for(int i=0;i<20;++i) {
    V g=V(std::sin(.71f*i),std::cos(.47f*i),std::sin(.29f*i+1)).normalized()*9.8025605f;
    V axis=V(std::sin(1.1f*i+1),std::cos(.91f*i+.2f),.8f).normalized();
    if(std::fabs(axis.dot(g.normalized()))>.85f)axis=g.normalized().unitOrthogonal();
    const V w=axis*(.8f+.04f*i),raw=S.inverse()*w+bias;
    fit->beginSegment(us,g,25,i%4==3);V end=g;
    for(int k=0;k<=200;++k) {
      us+=5000;fit->add(us,raw,25);
      if(k)end=Eigen::AngleAxisf(-w.norm()*.005f,axis)*end;
    }
    fit->finishSegment(end,25);us+=10000;
  }
  imu_cal::AccelCalibration<float> acc;acc.ok=acc.biasT.ok=true;acc.g=9.8025605f;
  imu_cal::GyroCalibration<float> gyro;gyro.ok=gyro.biasT.ok=true;gyro.biasT.b0=bias;
  int64_t start=esp_timer_get_time();const bool fitted=fit->fit(acc,gyro,[](void*){vTaskDelay(1);});
  Serial.printf("[BENCH] gyro_fit_ok=%d us=%lld buffer=%u samples=%d matrix_error=%.8g validation=%.8g stack_free=%uB\n",
    int(fitted),esp_timer_get_time()-start,unsigned(sizeof(Fit)),fit->samples(),
    double((gyro.S-S).norm()),double(gyro.rotation_validation_rms),unsigned(uxTaskGetStackHighWaterMark(nullptr)));
  fit->~Fit();heap_caps_free(memory);
  imu_cal::StationaryNoiseCapture<> noise;noise.begin(0);uint32_t worst=0;int64_t total=0;
  for(int i=0;i<=480001;++i) {
    V a(0,0,-9.8025605f),w;
    for(int k=0;k<3;++k){a[k]+=.02828427f*uniformNoise();w[k]=.002828427f*uniformNoise();}
    start=esp_timer_get_time();noise.add(uint32_t(i)*5000,a,w,25);
    const uint32_t dt=uint32_t(esp_timer_get_time()-start);total+=dt;worst=std::max(worst,dt);
    if(i%2048==0)vTaskDelay(1);
  }
  imu_cal::SensorNoiseModel report;const bool noise_ok=noise.report(report);
  Serial.printf("[BENCH] noise_ok=%d state=%uB average_us=%.6f worst_us=%u white_gyro=%.8g stack_free=%uB\n",
    int(noise_ok),unsigned(sizeof(noise)),double(total)/480002.,worst,double(report.white_density[3]),
    unsigned(uxTaskGetStackHighWaterMark(nullptr)));
  using Gate=ocean_imu::MagneticInnovationGate<float>;
  for(int mode=0;mode<3;++mode) {
    Gate gate;total=0;worst=0;
    const V r=mode==0?V(.1f,0,0):mode==1?V(std::sqrt(15.f),0,0):V(20,0,0);
    for(int i=0;i<1000;++i) {
      start=esp_timer_get_time();gate.evaluate(r,M::Zero(),M::Identity());
      const uint32_t dt=uint32_t(esp_timer_get_time()-start);total+=dt;worst=std::max(worst,dt);
    }
    Serial.printf("[BENCH] mag_mode=%d average_us=%.6f worst_us=%u state=%uB stack_free=%uB\n",
      mode,double(total)/1000.,worst,unsigned(sizeof(gate)),unsigned(uxTaskGetStackHighWaterMark(nullptr)));
  }
  Serial.println("[BENCH] done; synthetic timings do not qualify real sensor noise or calibration accuracy");
  vTaskDelete(nullptr);
}
void setup() {
  Serial.begin(115200);delay(1500);auto cfg=M5.config();M5.begin(cfg);
  if(xTaskCreatePinnedToCore(benchmark,"cal-bench",32768,nullptr,1,nullptr,0)!=pdPASS)
    Serial.println("[BENCH] worker allocation failed");
}
void loop(){delay(1000);}
