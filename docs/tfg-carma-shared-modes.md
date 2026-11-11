# Experimental shared-mode CARMA on the two-frame group

`KalmanTfgCarmaSharedModes<float>` combines the existing six-state-per-axis
minimal shared-mode CARMA(4,1) prior with the existing right-invariant TFG core.
It has 27 error states on `TwoFrameGroup<float,6,2>`. It is experimental; the
shipping OU-III estimator and its default selection are unchanged.

```cpp
#include "kalman_carma_shared_modes/SeaStateFusionFilter_TfgCarmaSharedModes.h"
using Fusion = ocean_imu::carma_shared_modes::SeaStateFusion_TfgCarmaSharedModes;
Fusion fusion;
Fusion::Config cfg;
// Set cfg sensor noise, gravity and magnetic service as for the CARMA wrapper.
fusion.begin(cfg); // Evidence excitation, estimated coupling and confidence
fusion.update(dt, gyro_body_ned, specific_force_body_ned, temperature_c);
fusion.updateMag(magnetometer_body_ned);
```

The state is `(R, [p,v,a,l,w1,w2], [bg,ba])`. Every motion column is a world
vector; both biases are body vectors. There is no artificial `S=0` update.
Six real realization states retain the CARMA(4,1) marginal law and complex
cross-axis phase relationships; they do not mean CARMA(6,3).

## Shared implementation

- `Kalman3D_Wave_TFG_CARMA3<T,NC,Motion>` retains its original two-parameter
  interface and independent model when `Motion=void`. The motion-policy path
  shares prediction, gyro/bias noise correlations, sensor updates, Joseph
  covariance correction, group injection and full covariance-reset transport.
- `CarmaSharedModesMotion.h` owns the existing minimal realization,
  input-normal factorization, stationary covariances, paired real/imaginary
  discretization, coupling estimator and model-coordinate maps. Its public
  reference-discretization accessors let TFG use the same algorithms.
- `SharedModesEstimator.h` is the existing confidence facade extracted from
  `KalmanCarmaSharedModes.h`. Both geometries use it. The old public name is
  an alias of the same QMEKF core plus this facade. A changed confidence
  rebuilds the excitation at the next step without waiting for a model
  commit, slow declines accumulate into covariance shrinks, attitude
  re-initialization restarts the coupling evidence, and copies keep the
  configuration.
- `kalman_common/ConditionalCovariance.h` contains the existing conditional Schur fallback,
  parameterized by index sets. TFG uses column-major world-vector indices;
  QMEKF retains axis-major indices. Both preserve the common-state-explained
  cross covariance and shrink only its conditional remainder.
- `SharedModesFusion<Core>` uses the existing CARMA schedule and startup.
  Its original QMEKF public name remains available. No second tuner or
  magnetic-acquisition policy is introduced.

The excitation keeps the vertical prior and the measured horizontal-on-vertical
regression `b`. When its applied coherence `b^H H^-1 b / sigma_z^2` would exceed
`SharedModesConfig::cap_coherence` (0.8), the horizontal prior is enlarged
along the coupling direction instead of shrinking `b`:
`H <- H + ((lambda - 1) / q) b b^H` with `q = b^H H^-1 b` and
`lambda = q / (0.8 sigma_z^2)`. The cross-wave prior is unchanged, and each
horizontal axis rescales its latent colouring by the square root of its
variance gain, exactly like an amplitude commit; `CouplingStatus` reports
those gains. The keel-frame horizontal correlation enters the baseline `H`.
Enlarging the whole block
instead inflates the cross-wave prior and costs horizontal accuracy where the
cap binds. `cap_coherence = 0` restores the earlier
cross-term shrinkage, which understates the regression and leaves only 0.1% of
the along-wave horizontal variance independent of heave when horizontal motion
exceeds the prior (a weak heave RAO). Where the
cap does not bind, nothing changes.

The motion process contribution between world columns is
`QR(i,j)*Re(Sigma) + QI(i,j)*Im(Sigma)`. TFG's gyro and gyro-bias cross terms
are added separately and retained. A model commit maps covariance through
physical additive coordinates before and after the per-axis latent maps.
This includes the attitude shear that a direct invariant row scaling misses.
No covariance-sized scratch matrix is added to the real-time stack, and all
sample-path matrices are fixed size.

Turning coupling off preserves the observable stationary CARMA process and
its lag covariance. Finite group corrections also transport the unobserved
coloring covariance, so an augmented nonlinear Gaussian TFG need not be
bitwise identical to its smaller independent realization after corrections.
The mathematical test reports that finite-reset sensitivity separately from
its strict process-law identities; the full-wrapper off comparison uses the
existing TFG/CARMA 1% equivalence criterion. No existing gate is relaxed.

## Reproduction

Use the repository's usual Eigen fallback, or pass `EIGEN_DIR=/path/to/eigen`.
The pinned dataset is the release used by the main workflow.

```sh
make all

# Reanalyse the committed compact scores without regenerating simulations.
python3 tools/tfg_shared_modes_study.py analyse
bash plots/kalman_tfg_shared_modes/draw_plots.sh
(cd doc/kalman_tfg_shared_modes && latexmk -lualatex -halt-on-error kalman_tfg_shared_modes.tex)

# A fresh directory avoids mixing the published executable/protocol hashes
# with a different compiler or later source revision.
TFG_STUDY_RESULTS=/tmp/tfg-shared-reproduction
# Primary 2x2 study: 24 seeds, 27 conditions, 1200 s at 200 Hz.
python3 tools/tfg_shared_modes_study.py run --jobs 5 --results "$TFG_STUDY_RESULTS"
# Existing no-transport variant and new no-transport variant, paired.
python3 tools/tfg_shared_modes_reset.py --jobs 2 --results "$TFG_STUDY_RESULTS"
# Eight pinned vessel records; conditional sensor/initialization resampling.
python3 tools/tfg_shared_modes_dataset.py --seeds 8 --jobs 2 --results "$TFG_STUDY_RESULTS"
python3 tools/tfg_shared_modes_transients.py --jobs 2 --results "$TFG_STUDY_RESULTS"
# Reduced-heave diagnosis: default cap, legacy cap and cap-target brackets.
python3 tools/tfg_shared_modes_weak_heave.py --jobs 2 --results "$TFG_STUDY_RESULTS"
python3 tools/tfg_shared_modes_study.py analyse --results "$TFG_STUDY_RESULTS"
# Copy the measured host_timing/host_jitter/host_fixed, ESP size/stack,
# firmware_builds, mathematical_verification and refactor_regression records
# into the fresh directory after repeating the measurement commands below.
python3 tools/tfg_shared_modes_report.py --results "$TFG_STUDY_RESULTS"
(cd doc/kalman_tfg_shared_modes && latexmk -lualatex -halt-on-error kalman_tfg_shared_modes.tex)
```

CI builds this article in the `kalman_tfg_shared_modes` matrix entry. It runs
the same plotting script, LaTeX compilation and retry, PDF artifact upload,
and release publication as the other articles, using the committed study
evidence without downloading or replaying the wave dataset.

`run` rejects caches from a different source/executable/protocol. Do not mix
runs after modifying an estimator. Independent paired seeds, score windows,
case definitions, excitation policy and source/executable hashes are recorded
in `reports/results/tfg_shared_modes/protocol.json`. Oracle cross spectra
appear only in names containing `ext-true` or another explicit external prior;
they never enter an operational default.
The original registered runner and implementation are preserved at the
revision in `execution_provenance.json`. Later runner changes permit
score-only reanalysis, fresh output directories and deterministic replay CSV
ordering; article claims are calculated from the measured result tables. They
do not alter the preserved experiment scores.
The registered revision preserves the simulation sources. Subsequent constructor
lint repairs use the same model order and covariance expressions; repeated
numerical equivalence checks are recorded in `review_validation.json`.
The existing OU-III finite secant
diagnostic binds the entire `src` tree, including experimental modules, so
it was rerun rather than weakening its verifier. Both perturbation sizes
reproduce the original raw shipping bytes exactly. Its refreshed source
binding and round-off differences are recorded in `shipping_faithfulness.json`;
all theorem qualifications remain unchanged.

The report keeps compact per-run scores (gzip), paired intervals, factorial
interactions, failures and one-second full-covariance audits. Raw 200 Hz
streams and restart caches are regenerated, not committed. Drawing families
have a 2000-point extrema-preserving budget, with the same additional vector
export decimation for SVG and PGF. Scores use the full-rate stream before
any drawing decimation.

The executable registry is shared with the existing CARMA benchmark:

| Label | Meaning |
|---|---|
| `A41@C` | Independent-axis QMEKF, Evidence excitation |
| `T41@C` | Independent-axis TFG, Evidence excitation |
| `CSM` | Existing shared QMEKF, estimated coupling and confidence |
| `TCSM` | New shared TFG, estimated coupling and confidence |
| `TCSM:off` | Same six-state realization with coupling disabled |
| `TCSM/real` | Estimated coupling with imaginary excitation removed |
| `TCSM/noreset` | Research ablation without reset covariance transport |
| `TCSM:noconf/ext-true` | Diagnostic truth-derived coupling reference |
| `TCSM/ext-flip` | Deliberately reversed coupling, confidence enabled |
| `TCSM:noconf/ext-flip` | Deliberately reversed coupling, unchecked |
| `TCSM/cap0`, `CSM/cap0` | Legacy PSD cap (cross term scaled down) |
| `TCSM/cap0.7`, `TCSM/cap0.9` | Cap-target brackets around the default 0.8 |
| `.../sf2` | Diagnostic horizontal prior factor 2 (default 1) |

The exploratory reduced-heave dose response reuses the same benchmark and the
first 12 primary seeds:

```sh
python3 tools/tfg_shared_modes_heave_sweep.py run --jobs 3
python3 tools/tfg_shared_modes_heave_sweep.py analyse   # committed scores only
```

It scales the heave RAO by `k` with horizontal motion unchanged and runs both
cap policies and both independent filters on the same data. It records the
applied multiple coherence of the known prior, and repeats the legacy-cap
cells near the boundary with the horizontal prior factor raised to 2. Its
repeated `k=1` and `k=0.15` cells must reproduce the primary and diagnosis
scores exactly.

The pinned replay is also shared source: `kalman_tfg_shared_modes-sim` compiles
`kalman_carma_shared_modes-sim.cpp` with the TFG core policy. Its results use a
distinct family label and do not reuse or modify the QMEKF's fitted regression
bars. The existing QMEKF replay retains its gating and public environment
variables. `W3D_CSM_VARIANT=CSM_OFF` selects the diagnostic uncoupled replay
for either binary. Independent-axis TFG can also be replayed with
`W3D_CARMA3_VARIANT=T41 tests/kalman_carma/carma3-sim`.

## Regression and mathematical verification

The baseline is main at `1fe4711b489ad92f3eef8ab560da904af741df9d`.
`tools/tfg_shared_modes_regression.py --baseline <pristine-checkout> --eigen
<Eigen-root>` compiles the same deterministic wrapper trace against both
header trees. It covers startup, changing excitation, bias handling, magnetic
updates and an outage; timings are excluded. Its results and the 27-state
independent-reference checks are committed beside the study.

The new mathematical target is `tfg_shared_modes-test`. It checks nonlinear
finite differences, independent continuous Lyapunov integration, full group
reset, per-axis adaptive coordinate shear, conditional PSD fallback, yaw
rotation equivariance, the outage yaw gauge, strict observable process-law
identities and float variable-interval operation under Eigen's allocation
guard. All existing mathematical and replay gates retain their thresholds.

## Computational measurements

```sh
make -C tests/kalman_carma_shared_modes kalman_carma_shared_modes-timing
./tests/kalman_carma_shared_modes/kalman_carma_shared_modes-timing --seconds 60
./tests/kalman_carma_shared_modes/kalman_carma_shared_modes-timing --seconds 60 --jitter
./tests/kalman_carma_shared_modes/kalman_carma_shared_modes-timing --seconds 60 --no-commits
EIGEN_DIR=/path/to/eigen bash tests/kalman_carma_shared_modes/arduino/csm_size_probe/build.sh /path/to/xtensa/bin /tmp/size-probe
python3 tools/tfg_shared_modes_firmware.py --build-root /tmp/tfg-shared-firmware
```

The timing driver runs four repetitions, rotating each core through all four
execution positions and recording `rep`/`order`. The committed `host_order_*`
files provide this sensitivity check separately from the original publication
timings; `review_validation.json` summarizes both.

The Xtensa probe covers all four factorial cores. Its stack-chain result is
a compiler-derived lower bound, not an observed peak stack. Code/literal and
read-only sizes are object-level measurements; full linked Arduino firmware
sizes are recorded separately. Host execution time does not establish ESP32
real-time deadlines. No new sample-path heap allocation is introduced; the
existing wrapper constructs its estimator on the heap at `begin()`.

The full AtomS3R application is the sketch
`sensors/full_marine_ins/atomS3R_ins_tfg_carma_shared_modes`: the TFG-CARMA3
sketch with the shared-mode TFG filter as its inner filter. The matched
firmware-size comparison (`tools/tfg_shared_modes_firmware.py`) instead
compiles the QMEKF shared-mode sketch with `-DOCEAN_IMU_TFG_SHARED_MODES=1` in
`compiler.cpp.extra_flags`, so that both builds share one application; that
sketch's default remains its QMEKF. The optional `tests/kalman_carma_shared_modes/arduino/tfg_shared_runtime_probe`
measures one selected core at 200 Hz for 60 s, including 0.1 s model commits
and 100 Hz magnetic updates. `PROBE_VARIANT=0..3` selects the four cells. It
prints mean/99th/max update time, missed 5 ms deadlines, minimum free heap,
and the FreeRTOS task's raw stack high-water value. It uses measured sample
intervals and needs no physical IMU. These measurements cover the core, not
the complete UI/sensor application. Target timing and stack high-water
validation still require an actual ESP32-S3 board; no physical board was
connected for this study.

The standalone article is
`doc/kalman_tfg_shared_modes/kalman_tfg_shared_modes.tex`. Its conclusions
use measured paired results, consistency and computational cost, including
small, unfavorable or insignificant differences.
