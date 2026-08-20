#!/usr/bin/env python3
"""Keep the deployed parameter-commit cadence out of the EMA-horizon experiment.

The experiment changes only averaging horizons.  The OU-III filter already
consumes every IMU/variance sample before this cadence gate; the gate merely
limits when the smoothed parameter candidate is activated.  Preserve it so the
sweep has no cadence confound.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

text = PATH.read_text()
old = """        // Every valid physical sample contributes to every EMA.  The candidate
        // is staged immediately and becomes active at the *next* IMU boundary,
        // preserving predictability without decimating measurements.
        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);
        online_tune_apply_pending_ = true;
    }
"""
new = """        // Every valid physical sample contributes to every EMA.  The deployed
        // activation cadence is deliberately kept separate from that sampling:
        // it throttles parameter commits, not physical measurements or EMA input.
        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);

        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            // y_k may change the smoothed candidate here, but the MEKF keeps
            // the schedule that was active before y_k arrived. Commit this
            // candidate at the beginning of updateTime(k+1).
            online_tune_apply_pending_ = true;
            last_adapt_time_sec_ = time_;
        }
    }
"""
if text.count(old) != 1:
    raise SystemExit(f"expected one experimental cadence block, found {text.count(old)}")
PATH.write_text(text.replace(old, new, 1))
print("PRESERVED deployed parameter-commit cadence; every EMA still consumes every sample")
