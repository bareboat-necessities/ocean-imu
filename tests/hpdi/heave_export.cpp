// heave_export -- per-sample export of one shipped estimator's vertical channel.
//
// Built once per estimator family (tests/hpdi/Makefile) from two units:
//   heave_export_shipped.cpp  the family's simulator source, included
//                             unchanged with main() renamed, plus accessors;
//   heave_export.cpp          this file: argument handling, hooks, writer.
// Nothing edits estimator code or the simulator sources; the replay evidence
// closure is untouched.
//
// Usage: heave_export-<family> --out EXPORT.csv --input RECORD.csv [simulator args]
// Every other argument and every W3D_* environment variable goes to the
// simulator unchanged (W3D_IMU_SEED, W3D_INIT_SEED, W3D_VALIDATION_WINDOW_SEC,
// W3D_TUNING_MODE, ...). W3D_VALIDATION_WINDOW_SEC must be set; set
// W3D_WRITE_TIMESERIES=0, as the studies do.
//
// OU-III, OU-II and TFG run the simulator's own main(). Two link-time hooks
// (-Wl,--wrap, see the Makefile) observe it without changing its code:
//   * the W3dSimulationRunner constructor is handed a Tap that forwards every
//     call to the shipped adapter and, for OU-III, reads that filter's front
//     end after each update;
//   * print_validation_metrics(), which the simulator calls with its run
//     result, keeps a copy of that result.
// The study compares this binary's VALIDATION_METRICS lines with the shipped
// simulator's byte for byte.
//
// The PII and TVG-NLO simulators score only their own deterministic records
// and print no VALIDATION_METRICS. Their builds call the simulator's own
// process function on the given record (seeds from the environment, as there)
// and print VALIDATION_METRICS with the shared print_validation_metrics().
//
// Output columns: z_ref, err_<F>, zhat_<F>; the OU-III build adds a_up and
// Tz_hat. One row per record sample, in record order (the study takes t from
// the record). err_<F> is the harness's own per-sample float error (estimate
// minus reference): exactly what its %Hs score accumulates. zhat_<F> is
// err_<F> + z_ref in double. a_up is the up-positive Mahony-proxy vertical
// acceleration (gravity removed, after the vibration guard) read from the
// OU-III filter's own front end after each sample; Tz_hat is that front end's
// canonical log-period estimate, NaN until its startup-usable gate clears.
// z_ref is the record's disp_z, which the harness scores as up-positive.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#define EIGEN_NON_ARDUINO
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "util/W3dSimCommon.h"

#if defined(HEAVE_EXPORT_OU_III)
#define HEAVE_EXPORT_FAMILY "OU_III"
#elif defined(HEAVE_EXPORT_OU_II)
#define HEAVE_EXPORT_FAMILY "OU_II"
#elif defined(HEAVE_EXPORT_TFG)
#define HEAVE_EXPORT_FAMILY "TFG"
#elif defined(HEAVE_EXPORT_PII)
#define HEAVE_EXPORT_FAMILY "PII"
#define HEAVE_EXPORT_DIRECT
#elif defined(HEAVE_EXPORT_TVG_NLO)
#define HEAVE_EXPORT_FAMILY "TVG_NLO"
#define HEAVE_EXPORT_DIRECT
#else
#error "define one HEAVE_EXPORT_<family>"
#endif

// Defined in heave_export_shipped.cpp.
int heave_export_shipped_main(int argc, char* argv[]);
#if defined(HEAVE_EXPORT_OU_III)
float heave_export_a_up(const IW3dFusionAdapter& adapter);
float heave_export_tz_hat(const IW3dFusionAdapter& adapter);
#endif
#if defined(HEAVE_EXPORT_DIRECT)
std::optional<W3dSimulationRunResult> heave_export_run_direct(const std::string& path, bool with_noise);
#endif

namespace heave_export {

struct Capture {
    std::optional<W3dSimulationRunResult> result;
    std::vector<float> a_up;
    std::vector<float> tz_hat;
};
Capture capture;

[[noreturn]] void die(const std::string& message) {
    std::fprintf(stderr, "heave_export: %s\n", message.c_str());
    std::exit(2);
}

class Tap final : public IW3dFusionAdapter {
public:
    explicit Tap(IW3dFusionAdapter& shipped) : shipped_(shipped) {}

    void updateMag(const Vector3f& mag_body_ned) override { shipped_.updateMag(mag_body_ned); }

    void update(float dt, const Vector3f& gyr, const Vector3f& acc, float temperature_c) override {
        shipped_.update(dt, gyr, acc, temperature_c);
#if defined(HEAVE_EXPORT_OU_III)
        capture.a_up.push_back(heave_export_a_up(shipped_));
        capture.tz_hat.push_back(heave_export_tz_hat(shipped_));
#endif
    }

    FilterSnapshot snapshot() const override { return shipped_.snapshot(); }

private:
    IW3dFusionAdapter& shipped_;
};

}  // namespace heave_export

// Link-time hooks. The by-value runner arguments are not trivial for the
// purposes of calls, so the Itanium ABI passes them by address; the hooks
// forward those addresses unchanged.
static_assert(!std::is_trivially_destructible_v<W3dSimulationOptions>);
static_assert(!std::is_trivially_destructible_v<SimulationNoiseModels>);

#define HEAVE_EXPORT_RUNNER_CTOR \
    _ZN19W3dSimulationRunnerC1E20W3dSimulationOptions21SimulationNoiseModelsR17IW3dFusionAdapter
#define HEAVE_EXPORT_PRINT_METRICS _Z24print_validation_metricsRK22W3dSimulationRunResultffPKc
#define HEAVE_EXPORT_STR2(x) #x
#define HEAVE_EXPORT_STR(x) HEAVE_EXPORT_STR2(x)
#define HEAVE_EXPORT_CAT2(a, b) a##b
#define HEAVE_EXPORT_CAT(a, b) HEAVE_EXPORT_CAT2(a, b)

void heave_export_real_runner_ctor(W3dSimulationRunner* self, W3dSimulationOptions* options,
                                   SimulationNoiseModels* noise, IW3dFusionAdapter* adapter)
    asm(HEAVE_EXPORT_STR(HEAVE_EXPORT_CAT(__real_, HEAVE_EXPORT_RUNNER_CTOR)));
void heave_export_real_print_metrics(const W3dSimulationRunResult& result, float dt, float window,
                                     const char* family)
    asm(HEAVE_EXPORT_STR(HEAVE_EXPORT_CAT(__real_, HEAVE_EXPORT_PRINT_METRICS)));

extern "C" void HEAVE_EXPORT_CAT(__wrap_, HEAVE_EXPORT_RUNNER_CTOR)(W3dSimulationRunner* self,
                                                                   W3dSimulationOptions* options,
                                                                   SimulationNoiseModels* noise,
                                                                   IW3dFusionAdapter* adapter) {
    static std::optional<heave_export::Tap> tap;
    if (tap) heave_export::die("one record per process");
    tap.emplace(*adapter);
    heave_export_real_runner_ctor(self, options, noise, &*tap);
}

extern "C" void HEAVE_EXPORT_CAT(__wrap_, HEAVE_EXPORT_PRINT_METRICS)(const W3dSimulationRunResult& result,
                                                                     float dt, float window,
                                                                     const char* family) {
    heave_export::capture.result = result;
    heave_export_real_print_metrics(result, dt, window, family);
}

namespace heave_export {

int run_main(int argc, char** argv) {
    std::string input, output;
    bool with_noise = true;
    std::vector<char*> forwarded{argv[0]};
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--out" && i + 1 < argc) {
            output = argv[++i];
            continue;
        }
        if (arg == "--input" && i + 1 < argc) {
            if (!input.empty()) die("exactly one --input is supported");
            input = argv[i + 1];
        }
        if (arg == "--no-noise") with_noise = false;
        forwarded.push_back(argv[i]);
    }
    if (input.empty() || output.empty()) die("usage: --out EXPORT.csv --input RECORD.csv [simulator args]");
    float window_sec = 0.0f;
    if (const char* s = std::getenv("W3D_VALIDATION_WINDOW_SEC")) window_sec = static_cast<float>(std::atof(s));
    if (!(window_sec > 0.0f)) die("W3D_VALIDATION_WINDOW_SEC must be positive");

    int status = 0;
#if defined(HEAVE_EXPORT_DIRECT)
    (void)forwarded;
    std::optional<W3dSimulationRunResult> direct;
    try {
        direct = heave_export_run_direct(input, with_noise);
    } catch (const std::exception& e) {
        die(e.what());
    }
    if (direct) print_validation_metrics(*direct, 1.0f / 200.0f, window_sec, HEAVE_EXPORT_FAMILY);
#else
    (void)with_noise;
    forwarded.push_back(nullptr);
    status = heave_export_shipped_main(static_cast<int>(forwarded.size()) - 1, forwarded.data());
#endif
    if (!capture.result) die("record was not processed: " + input);
    const W3dSimulationRunResult& r = *capture.result;
    const size_t n = r.errs_z.size();
    if (r.ref_z.size() != n) die("reference and error lengths differ");
#if defined(HEAVE_EXPORT_OU_III)
    if (capture.a_up.size() != n || capture.tz_hat.size() != n) die("front-end tap is misaligned");
#endif

    FILE* fp = std::fopen(output.c_str(), "w");
    if (!fp) die("cannot write " + output);
    std::fprintf(fp, "z_ref,err_%s,zhat_%s", HEAVE_EXPORT_FAMILY, HEAVE_EXPORT_FAMILY);
#if defined(HEAVE_EXPORT_OU_III)
    std::fprintf(fp, ",a_up,Tz_hat");
#endif
    std::fprintf(fp, "\n");
    for (size_t k = 0; k < n; ++k) {
        const double z_hat = static_cast<double>(r.errs_z[k]) + static_cast<double>(r.ref_z[k]);
        std::fprintf(fp, "%.9g,%.9g,%.10g", static_cast<double>(r.ref_z[k]), static_cast<double>(r.errs_z[k]),
                     z_hat);
#if defined(HEAVE_EXPORT_OU_III)
        std::fprintf(fp, ",%.9g,%.9g", static_cast<double>(capture.a_up[k]),
                     static_cast<double>(capture.tz_hat[k]));
#endif
        std::fprintf(fp, "\n");
    }
    std::fclose(fp);
    std::printf("HEAVE_EXPORT family=%s samples=%zu hs=%.9g\n", HEAVE_EXPORT_FAMILY, n,
                static_cast<double>(r.wave_params.height));
    return status;
}

}  // namespace heave_export

int main(int argc, char** argv) { return heave_export::run_main(argc, argv); }
