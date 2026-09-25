# Wrapper equivalence traces

Bit-level characterization of the three sea-state orchestrators
(`SeaStateFusion_OU_II`, `SeaStateFusion_OU_III`, `SeaStateFusionFilter_TFG`
and the two inner OU filters), for refactors that must not change behavior.

`SyntheticMarineHistory.h` generates deterministic 200 Hz histories: wave-band
roll/pitch and orbital acceleration, gyro/accelerometer bias, a hard-iron
offset, 25 Hz magnetometer, an engine-vibration interval, a still-water
interval and an optional capsize.  The drivers read only public API and write
every sample's outputs as raw float bits:

* attitude (wrapper and MEKF), velocity, displacement, integral (OU-III/TFG),
  world acceleration, gyro and accelerometer bias, covariance diagonal and
  selected cross terms;
* tau / sigma_aw / regularizer targets and applied values, pseudo-update
  period, information ratio;
* startup stage, tuner readiness, live/north-lock/refine times, magnetic
  acceptance, hard-iron startup and continuous offsets;
* tracker frequencies, wave period, sigma band, direction axis, sign state,
  confidence and coherence;
* vibration-guard engagement, removed RMS, commanded Racc;
* displacement detrender output.

The inner-filter scenarios also flip runtime knobs (fixed tuning, law
switches, cadence toggles, congruent covariance sync, channel freezes, bounds,
low-wave noise weighting) and drive the discouraged pre-Live path through a
capsize to reach the Cold re-entry.  Each run prints a coverage line.

Usage, comparing a branch against `main`:

    git worktree add /tmp/base origin/main
    tools/wrapper_equivalence/build_and_trace.sh /tmp/base/src /tmp/traces-base
    tools/wrapper_equivalence/build_and_trace.sh src            /tmp/traces-new
    python3 tools/wrapper_equivalence/compare_traces.py /tmp/traces-base /tmp/traces-new

The default build uses `-ffp-contract=off`, so a trace depends only on the
order of floating-point operations: a behavior-preserving refactor must be
byte-identical.  With FMA contraction enabled (the deployed `-march=native`
flags, via `CXXFLAGS=...`) moving code between functions can change which
products the compiler fuses and therefore low-order bits, without changing a
single source operation; compare_traces.py reports the size of such drift.
