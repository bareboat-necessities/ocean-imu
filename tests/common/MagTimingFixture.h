#pragma once
// Synthetic qualified timing metadata for persistence tests; never device evidence.
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
inline void bindTimingFixture(atoms3r_ical::ImuCalBlobV5& b, uint32_t configuration=123) {
    using namespace atoms3r_ical;
    auto& r=b.mag_timing;
    r={}; r.binding=magTimingBinding(b); r.configuration=configuration;
    r.delay_us=40000; r.fit_width_s=.002f; r.residual_time_s=.025f;
    r.train_rms_rad=r.validation_rms_rad=.005f; r.improvement=.8f;
    r.train_pairs=r.validation_pairs=100; r.crc=magTimingCrc(r);
}
