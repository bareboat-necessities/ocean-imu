/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R Tilt-Compensated Compass + IMU Calibration Wizard (Mahony)
*/

#include <Arduino.h>
#include <M5Unified.h>

#define NO_BOSCH_API

// 1 = graphical compass by default, 0 = text UI by default
#ifndef COMPASS_UI_DEFAULT_GRAPHICS
#define COMPASS_UI_DEFAULT_GRAPHICS 1
#endif

// 0 = keep debug serial
// 1 = emit NMEA0183 (HDM + XDR + ROT) like pypilot
#ifndef COMPASS_SERIAL_NMEA
#define COMPASS_SERIAL_NMEA 1
#endif

#ifndef COMPASS_NMEA_TALKER
#define COMPASS_NMEA_TALKER "II"
#endif

#include "AtomS3R/AtomS3R_CompassAppBase.h"

using namespace atoms3r_compass;

class MahonyCompassApp : public CompassAppBase {
 public:
  MahonyCompassApp()
      : CompassAppBase(std::make_unique<MahonyBackend>(), MagGateConfig{12, 0.02f, 10, 20}, "Mahony") {}
};

static MahonyCompassApp g_app;

void setup() { g_app.begin(); }
void loop() { g_app.tick(); }
