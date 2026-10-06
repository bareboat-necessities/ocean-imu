// Shipped half of heave_export (see heave_export.cpp).
//
// This translation unit is one simulator source, included unchanged with its
// main() renamed, plus a few small accessors. Keep it that way. GCC sizes its
// inlining budget by the whole unit and contracts floating-point expressions
// into FMAs according to those decisions, so substantial extra code here
// changes the replay in the last bits. The study compares this build's
// VALIDATION_METRICS lines with the shipped simulator's byte for byte.

#define main heave_export_shipped_main
#if defined(HEAVE_EXPORT_OU_III)
#include "../kalman_ou_iii/kalman_ou_iii-sim.cpp"
#elif defined(HEAVE_EXPORT_OU_II)
#include "../kalman_ou_ii/kalman_ou_ii-sim.cpp"
#elif defined(HEAVE_EXPORT_TFG)
#include "../kalman_tfg/kalman_tfg-sim.cpp"
#elif defined(HEAVE_EXPORT_PII)
#include "../pii_observer/pii_observer-adaptive.cpp"
#elif defined(HEAVE_EXPORT_TVG_NLO)
#include "../nlo/nlo-sim.cpp"
#else
#error "define one HEAVE_EXPORT_<family>"
#endif
#undef main

#if defined(HEAVE_EXPORT_OU_III)
// FusionAdapter_OU_III keeps its SeaStateFusion_OU_III private and the
// simulator source must not change. An explicit instantiation is exempt from
// access checking ([temp.spec.general]/6), which lets this file name that
// member without editing the class.
using HeaveExportOU3Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
template <typename Tag, typename Tag::type Member>
struct HeaveExportPrivateMember {
    friend typename Tag::type heave_export_member(Tag) { return Member; }
};
struct HeaveExportOU3FusionTag {
    using type = HeaveExportOU3Fusion FusionAdapter_OU_III::*;
    friend type heave_export_member(HeaveExportOU3FusionTag);
};
template struct HeaveExportPrivateMember<HeaveExportOU3FusionTag, &FusionAdapter_OU_III::fusion_>;

// The adapter the runner steps is the simulator's FusionAdapter_OU_III; these
// read the measurement-only front end of its filter.
static const auto& heave_export_filter(const IW3dFusionAdapter& adapter) {
    const auto& shipped = static_cast<const FusionAdapter_OU_III&>(adapter);
    return (shipped.*heave_export_member(HeaveExportOU3FusionTag{})).raw();
}
float heave_export_a_up(const IW3dFusionAdapter& adapter) {
    return heave_export_filter(adapter).getAccelVertical();
}
float heave_export_tz_hat(const IW3dFusionAdapter& adapter) {
    const auto& filter = heave_export_filter(adapter);
    return filter.wavePeriodUsable() ? filter.getWavePeriodSec() : NAN;
}
#endif

#if defined(HEAVE_EXPORT_PII)
std::optional<W3dSimulationRunResult> heave_export_run_direct(const std::string& path, bool with_noise) {
    add_noise = with_noise;
    return process_wave_file_for_adaptive_pii_mahony(path, 1.0f / 200.0f, true, add_noise, 20.0f);
}
#elif defined(HEAVE_EXPORT_TVG_NLO)
std::optional<W3dSimulationRunResult> heave_export_run_direct(const std::string& path, bool with_noise) {
    add_noise = with_noise;
    auto result = process_wave_file_for_tvg_nlo_nomag_nognss(path, 1.0f / 200.0f, add_noise);
    if (!result) return std::nullopt;
    return static_cast<const W3dSimulationRunResult&>(*result);
}
#endif
