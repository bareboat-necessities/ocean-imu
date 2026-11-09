// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Flash-only measurement harness, no physical IMU needed. One core at 200 Hz.
// PROBE_VARIANT: 0 independent QMEKF, 1 independent TFG, 2 shared QMEKF,
// 3 shared TFG. This measures the estimator, not the complete UI/sensor app.
#include <Arduino.h>
#include <algorithm>
#include <ArduinoEigenDense.h>
#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#ifndef PROBE_VARIANT
#define PROBE_VARIANT 3
#endif
#if PROBE_VARIANT == 0
using Core=ocean_imu::kalman::Kalman3D_Wave_CARMA3<float,4,false,4,false>;
#elif PROBE_VARIANT == 1
using Core=ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float,4>;
#elif PROBE_VARIANT == 2
using Core=ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>;
#else
using Core=ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>;
#endif
using V=Eigen::Vector3f;
constexpr unsigned count=12000;
uint32_t timings[count]; // Drawing/quantile storage is outside the task stack.
Core* core=nullptr;
TickType_t wake_tick;
uint32_t previous_start=0;
unsigned sample=0,deadline_misses=0;
uint64_t sum_us=0;
void setup(){
    Serial.begin(115200);delay(1500);
    core=new Core(V::Constant(.0294f),V::Constant(.000031f),V::Constant(7.7f));
    core->set_vertical_carma(1.f,.6f);core->set_horizontal_carma(1.f,.5f,.3f);
    core->reset_aw_covariance_to_stationary();
    core->initialize_from_acc_mag(V(0,0,-9.80665f),V(20,0,45));
    wake_tick=xTaskGetTickCount();
}
void loop(){
    vTaskDelayUntil(&wake_tick,pdMS_TO_TICKS(5));
    const uint32_t now=micros();const float dt=previous_start?float(now-previous_start)*1e-6f:.005f;
    previous_start=now;const float t=float(sample)*.005f;
    const V acc(.3f*std::sin(1.9f*t),.2f*std::cos(1.9f*t),-9.80665f+.4f*std::cos(1.9f*t));
    const uint32_t start=micros();
    if(sample%20==0){const float k=1.f+.1f*std::sin(.1f*t);core->set_vertical_carma(k,.6f);core->set_horizontal_carma(k,.5f,.3f);}
    core->time_update(V(.001f,-.002f,.0005f),dt);
    core->measurement_update_acc_only(acc,35.f);
    if(sample%2==0)core->measurement_update_mag_only(V(20,0,45));
    const uint32_t elapsed=micros()-start;timings[sample]=elapsed;sum_us+=elapsed;deadline_misses+=elapsed>5000;
    if(++sample==count){
        // Capture before quantile sorting/printing adds unrelated stack frames.
        const auto highwater=uxTaskGetStackHighWaterMark(nullptr);
        const uint32_t minheap=ESP.getMinFreeHeap();
        std::sort(timings,timings+count);
        Serial.printf("variant,states,core_bytes,mean_us,p99_us,max_us,missed_5ms,task_stack_highwater,min_free_heap\n");
        Serial.printf("%d,%d,%u,%.3f,%lu,%lu,%u,%lu,%lu\n",PROBE_VARIANT,Core::NX,unsigned(sizeof(Core)),
                      double(sum_us)/count,(unsigned long)timings[count*99/100],(unsigned long)timings[count-1],
                      deadline_misses,(unsigned long)highwater,(unsigned long)minheap);
        // The platform's high-water units follow its FreeRTOS port. Retain the
        // raw value and board/core version; do not infer bytes in host reports.
        for(;;)vTaskDelay(pdMS_TO_TICKS(1000));
    }
}
