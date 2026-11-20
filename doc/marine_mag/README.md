# Marine magnetic calibration study

Read the standalone [article](marine-mag.tex) and the [paired evidence](../../reports/results/marine_mag/). The baseline is f1424f58ef019994a28aef03efc3da0dc5d5ca98. Each navigation run processes a 1200 s, 200 Hz release record or continuous regenerated maneuver and scores its final 900 s.

### Heading improvement against the previous always-on filter

TFG shared-mode CARMA now uses a 60-second calibration slew. Each calibration step explains part of the change in the measured horizontal field bearing. The filter applies 20% of that rotation about world vertical before the existing magnetic update. The normal simulator uses this response automatically. `W3D_MAG_CALIBRATION=native` reproduces the earlier 45-second response, and `production` selects the current filter configuration. The saved-calibration sketch still controls whether continuous learning is enabled.

| Panel | Paired records | Previous yaw RMS | New yaw RMS | Paired change (95% CI) |
|---|---:|---:|---:|---:|
| Released records | 16 | 0.731° | 0.711° | -0.020°  |
| 8 fresh sensor draws | 48 | 0.959° | 0.931° | -0.028° [-0.062, +0.007] |
| 32 fresh sensor draws | 192 | 0.789° | 0.767° | -0.022° [-0.036, -0.009] |

On the 16 released records, following-sea yaw RMS went from 0.760° to 0.747° and head-sea RMS from 0.702° to 0.676°. All original navigation gates passed on both default and interval-rate profiles. The only yaw increase was 0.0012° in the 8.5 m PM–Stokes head sea. The focused replay also checks a mean yaw ceiling independently for each direction.

The response was selected using the released records. Independent draws vary sensor noise and installation/initialization together. Each draw covers PM–Stokes at 0.27, 1.5 and 8.5 m in both directions. The simulator applies its per-record limits to the default seed. We evaluate fresh draws by paired yaw RMS. Intervals use paired draw means. The new response improved yaw in 7 of 8 validation draws and 22 of 32 subsequent replication draws. The earlier fixed-time switch that gave 0.626° on released records gave 1.109° on four fresh draws, compared with 1.040° for the previous learner.

~~~sh
make -C tests/kalman_tfg_shared_modes kalman_tfg_shared_modes-sim
for panel in released validation replication; do
  python3 tools/marine_mag_yaw_response.py "$panel" --out "/tmp/yaw-$panel.csv"
done
~~~

[Released records](../../reports/results/marine_mag/yaw-response-released.csv) · [Eight-draw validation](../../reports/results/marine_mag/yaw-response-validation.csv) · [32-draw replication](../../reports/results/marine_mag/yaw-response-replication.csv) · [Timing comparison](../../reports/results/marine_mag/yaw-response-timing-ablation.csv) · [Interval-rate gates](../../reports/results/marine_mag/yaw-response-interval-rate.csv)

### Installation screening

On eight fresh draws, the installation screen reduced mean sustained yaw RMS on clean noisy inputs from 0.765° to 0.469° against always-on learning. Distorted inputs moved from 0.898° to 0.907°. Eight further draws in 0.27 m waves reproduced the clean benefit (0.624° to 0.392°), while distorted inputs moved from 0.778° to 0.801°. The largest delayed-release penalty was 1.182°. The production implementation is available as an opt-in. [Primary paired effects](../../reports/results/marine_mag/installation-screen-paired.csv) · [Weak-wave confirmation](../../reports/results/marine_mag/installation-screen-weak-paired.csv)

The committed research scores use Debian 13 glibc 2.41 and the repository's replay compiler. The focused CI uses that userspace because different libm implementations can shift the final digits. Run from the repository root with the replay compiler, Eigen, GNU patch, NumPy, SciPy, pandas, Matplotlib and TeX Live.

~~~sh
for side in following head; do
  OCEAN_IMU_SIM_CASE="$side" python3 tools/sim_dataset.py --dest ".sim-data/cases/$side"
done
make -C tests/kalman_tfg_shared_modes/magnetic_study all check
export OPENBLAS_NUM_THREADS=1
python3 tools/marine_mag_information.py
for stage in diagnosis development validation stress; do
  python3 tools/marine_mag_study.py "$stage" --cache /tmp/ocean-magnetic-results --jobs 4
done
make -C tests/kalman_tfg_shared_modes/magnetic_study build/magnetic-quarantine-study
python3 tools/marine_mag_quarantine.py
make -C tests/kalman_tfg_shared_modes/magnetic_study build/magnetic-heading-study
python3 tools/marine_mag_yaw.py
python3 tools/marine_mag_yaw.py --gyro-bias
python3 tools/marine_mag_yaw.py --replicate-gyro-bias
python3 tools/marine_mag_yaw.py --gyro-stress
python3 tools/marine_mag_attribution.py
python3 tools/marine_mag_attribution.py --weak
python3 tools/marine_mag_attribution.py --stress
python3 tools/marine_mag_reference.py
python3 tools/marine_mag_transients.py
python3 tools/marine_mag_report.py
python3 tests/kalman_tfg_shared_modes/magnetic_study/check_evidence.py
(cd doc/marine_mag && latexmk -pdf marine-mag.tex)
~~~

The experiment builds a private include tree under the ignored build directory. Its instrumentation patch applies with zero fuzz to copies of baseline sources. It gives a named research adapter access to the posterior, reuses packet prediction during quarantine, and substitutes controlled sensor fixtures. The stock calibration equations, wrapper, reference slew and TFG update supply all arms. The isolated variants are diagnostic. The production integration below provides the screened installation setting. Use the magnetic-study executable's positional arguments (CSV, variant, scenario, dip in degrees). Arm -1 uses the sketch calibration switch. V0 through V8 are defined in the article. V9 is the supplemental V6-plus-covariance development arm. The separate magnetic-quarantine-study binary applies V8 packet rejection with the native measurement covariance. Its 32 post-hoc paired results are in quarantine-only.csv.

Each stage publishes a SHA-256 input manifest. Weak-steering and larger-turn tests replace the release track and rebuild the IMU. A companion track file is read only by the displacement scorer. W3D expands the joint seed into independent sensor and initialization streams. The stress manifest also hashes the four maneuver scoring tracks. Resumable checkpoints key the physical and scoring inputs, executable, scoring script and experimental condition.

Navigation scores use the full sample rate. Selected figures retain uniform five-second summaries at ten-second spacing. Calibration RMS covers five-second epochs from 305 through 1195 s. Packet counters and information are recorded through 1195 s. `accepted` counts native EKF updates, `rejected` counts external quarantine events, and `missed` counts forwarded post-initialization packets for which the native accepted count stayed constant. An unapplied correction has first_applied = -1. Confidence intervals use one paired effect per independent draw. Rows from the same draw share initialization and sensor noise across directions.

The optional magnetic-reference-study executable imposes an external map at the core-update boundary or injects an error at the live handoff. The primary calibration variants cannot select these oracle controls.

For an ESP32-S3 cross-compile, supply the bin directory of Espressif GCC esp-14.2.0_20241119 to the command below. It emits object sizes, static stack diagnostics and a separate host microbenchmark. Board runtime timing and stack high-water require a physical device.

~~~sh
bash tools/marine_mag_resources.sh /path/to/xtensa-esp-elf/bin /tmp/marine-mag-resources
~~~

For a separate sanitizer build, run the following after fetching both datasets. On a host that disallows LeakSanitizer under ptrace, set `ASAN_OPTIONS=detect_leaks=0`.

~~~sh
study=tests/kalman_tfg_shared_modes/magnetic_study
sanout=$(mktemp -d)
make -j1 -C "$study" OUT="$sanout" CXX=clang++ \
  CXXFLAGS="-std=c++20 -O1 -g0 -DEIGEN_UNROLLING_LIMIT=0 -fsanitize=address,undefined -fno-omit-frame-pointer" all check
UBSAN_OPTIONS=halt_on_error=1 "$sanout/magnetic-study" .sim-data/cases/following/wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv 8 burst 66.5
~~~

The sustained-heading extension compares post-refinement magnetic standard deviation factors of 0.5, 1, 2 and 4. The development screen requires lower average yaw RMS and absolute signed yaw, mean displacement growth at most 2% with a 10% worst-case cap, and mean roll/pitch RMS growth at most 0.01 degrees with a 0.1-degree worst-case cap. If every alternative fails, the lowest-yaw development candidate is tested for falsification. The fresh 32-cell panel uses four new seeds and uninterrupted records. yaw-weight-paired.csv reports paired draw-level intervals, counts and worst changes.

The gyro-bias experiment uses the same development screen with the native TFG magnetic weight. It scales the external learner's initial gyro-bias standard deviation and gyro-bias random-walk amplitude together by 0, 0.25, 1 or 4. The zero factor freezes that nuisance state at its initial zero. The chosen alternative and the native learner are evaluated on four further joint draws (547, 1151, 2371 and 4099). The three gyro-bias CSVs contain the scored arms and their paired comparisons.

The independent confirmation panel fixes the native and fourfold gyro-bias factors before evaluation, keeping the same eight uninterrupted physical conditions. Twelve further independent joint draws test transfer across sensor and initialization errors. Evidence for a sustained yaw improvement requires the upper endpoints of the paired 95% intervals for both total yaw RMS and absolute within-record bias to fall below zero on this panel, together with the development roll, pitch and displacement limits.

The installation screen crosses three magnetic inputs (ideal, white noise only, and white noise plus installation distortion) with three settings (sketch opt-out, native continuous calibration, and gated application V11). It runs the same eight release records under each of eight previously unused joint draws. V11 estimates distortion at the most recent posterior every five seconds. Three successive scores above 16.27 release the native application slew. The input generator consumes the same random draws in every arm. The candidate uses only the calibration posterior and its covariance. The two installation-screen CSVs contain 576 per-run scores and paired draw-level intervals. A V11 release time of -1 means that application stayed gated for the entire record.

The independent weak-motion confirmation fixes V11 and the 0.27 m wave record, changes dip between 45 and 76 degrees in both directions, and repeats all nine input/policy combinations on eight new draws. The separate stress run repeats the original 14 stress fixtures for V0 and V11.

### Production integration

The proxy-startup filters and the standalone TFG wrapper accept `mag_screen_installation` alongside `mag_continuous_hard_iron`. Enable both with the gyro-aided calibration to hold application until the posterior distortion passes three five-second checks. The shared-mode TFG–CARMA AtomS3R sketch enables both by default. Build it with `-DOCEAN_IMU_SCREEN_MAG_INSTALLATION=0` to disable online correction. The six other saved-calibration sketches that support screening opt in with `-DOCEAN_IMU_SCREEN_MAG_INSTALLATION=1`. The ordinary navigation filters use elapsed seconds to support varying packet rates. The original V11 evidence uses 125 received packets per check at 25 Hz. Both paths use the same posterior distortion calculation.

Build the normal simulator and run `python3 tools/marine_mag_filter_replays.py tfg_shared --jobs 4 --out /tmp/tfg_shared.csv` to compare disabled, native and screened continuous calibration on all 16 released records. The driver records individual scores and the original simulator exit status for every arm. The production comparison workflow runs the eight filter families and publishes their CSVs.

### Production replay results

Each value is the mean of 16 per-record yaw RMS values. Both directions, both spectra and four sea states contribute equally. Each score uses the final 900 seconds. Across the 128 always-on records, all gate results matched independent pre-integration reruns. The largest difference in yaw RMS was 0.000076 degrees. 112 records matched all six reported metrics exactly.

| Filter | Off (°) | Always on (°) | Screened (°) | Gate passes out of 16 (off/on/screened) |
| --- | ---: | ---: | ---: | ---: |
| TFG shared-mode CARMA | 2.732 | 0.731 | 0.786 | 0/16/12 |
| QMEKF shared-mode CARMA | 2.735 | 0.742 | 0.798 | 0/16/12 |
| Tangent shared-mode CARMA | 2.738 | 0.736 | 0.791 | 0/16/12 |
| TFG CARMA | 2.738 | 0.735 | 0.790 | 0/16/12 |
| QMEKF CARMA | 2.741 | 0.747 | 0.802 | 0/16/12 |
| OU-II | 2.830 | 0.816 | 0.865 | 0/16/12 |
| OU-III | 2.776 | 0.759 | 0.813 | 0/16/12 |
| TFG | 2.576 | 0.684 | 0.723 | 0/16/10 |

Screening lowered yaw RMS against correction off in all 128 comparisons. The complete quality gates passed 0/128 with correction off, 128/128 always on and 94/128 screened. Each of the 34 screened failures exceeded its existing gyro-bias limit.

For TFG shared-mode CARMA, screening reduced mean yaw RMS by 71% against correction off. With a distorted installation, always-on correction yielded the lowest average and passed all 16 quality gates. On all 12 records above 0.27 m, screened and always-on yaw RMS agreed at the printed precision. The four 0.27 m records exceeded the gyro-bias gate under screening. Mean 3-D displacement RMS was 0.15068 m off, 0.14754 m always on and 0.14754 m screened. On clean noisy installations the independent study above favors screening.

[Hosted comparison run](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44665619717) · [All 384 per-record scores, gate results and reasons](../../reports/results/marine_mag/production-filter-replays.csv)
