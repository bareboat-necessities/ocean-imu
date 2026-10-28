/*
  Link-only driver for the size probe (callgraph.sh): begin(), update() and
  updateMag() through the probe's C entry points, so every symbol the
  firmware sample path needs must resolve when the probe object is linked
  against the Xtensa toolchain's newlib/libstdc++.  Never run.
*/

// Built only by callgraph.sh; host compilers see an empty translation unit
// (as carma3_size_probe.cpp).
#if defined(__XTENSA__)

extern "C" void* carma3_probe_make();
extern "C" float carma3_probe_step(void* f, float dt, const float* gyr, const float* acc, const float* mag);

volatile float carma3_probe_sink;

int main() {
    void* f = carma3_probe_make();
    const float gyr[3] = {0.001f, -0.002f, 0.0005f};
    const float acc[3] = {0.1f, -0.05f, -9.8f};
    const float mag[3] = {20.0f, 1.0f, 40.0f};
    for (int k = 0; k < 1000; ++k) carma3_probe_sink = carma3_probe_step(f, 0.005f, gyr, acc, (k % 8) ? nullptr : mag);
    return 0;
}

#endif  // __XTENSA__
