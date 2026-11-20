# Marine magnetic calibration study

Read the standalone [article](marine-mag.tex) and the [paired evidence](../../reports/results/marine_mag/). The baseline is f1424f58ef019994a28aef03efc3da0dc5d5ca98. Each navigation run processes a 1200 s, 200 Hz release record or continuous regenerated maneuver and scores its final 900 s.

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
python3 tools/marine_mag_reference.py
python3 tools/marine_mag_transients.py
python3 tools/marine_mag_report.py
python3 tests/kalman_tfg_shared_modes/magnetic_study/check_evidence.py
(cd doc/marine_mag && latexmk -pdf marine-mag.tex)
~~~

The experiment builds a private include tree under the ignored build directory. Its instrumentation patch applies with zero fuzz to copies of baseline sources. It gives a named research adapter access to the posterior, reuses packet prediction during quarantine, and substitutes controlled sensor fixtures. The stock calibration equations, wrapper, reference slew and TFG update supply all arms. The calibration research adds no firmware setting. Use the magnetic-study executable's positional arguments (CSV, variant, scenario, dip in degrees). Arm -1 uses the sketch calibration switch. V0 through V8 are defined in the article. V9 is the supplemental V6-plus-covariance development arm.

Each stage publishes a SHA-256 input manifest. Weak-steering and larger-turn tests replace the release track and rebuild the IMU. A companion track file is read only by the displacement scorer. W3D expands the joint seed into independent sensor and initialization streams. Resumable checkpoints key the scored fixture, physical input and executable. Use a fresh cache after changing scoring code.

Navigation scores use the full sample rate. Selected figures retain uniform five-second summaries at ten-second spacing. Calibration RMS covers five-second epochs from 305 through 1195 s. Packet counters and information are recorded through 1195 s. An unapplied correction has first_applied = -1. Confidence intervals use one paired effect per independent draw. Rows from the same draw share initialization and sensor noise across directions.

The optional magnetic-reference-study executable imposes an external map at the core-update boundary or injects an error at the live handoff. The primary calibration variants cannot select these oracle controls.

For an ESP32-S3 cross-compile, supply the bin directory of Espressif GCC esp-14.2.0_20241119 to the command below. It emits object sizes, static stack diagnostics and a separate host microbenchmark. Board runtime timing and stack high-water require a physical device.

~~~sh
bash tools/marine_mag_resources.sh /path/to/xtensa-esp-elf/bin /tmp/marine-mag-resources
~~~
