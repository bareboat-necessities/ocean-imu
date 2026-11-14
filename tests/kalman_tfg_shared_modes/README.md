# TFG shared-mode CARMA deterministic regression

`kalman_tfg_shared_modes-sim.cpp` instantiates the 27-state TFG estimator with
`tests/common/SharedModesReplay.h`. The QMEKF uses that same adapter and replay
protocol. The common CARMA runner supplies the sensor streams, startup and metrics.

```sh
make -C tests/kalman_tfg_shared_modes all
# Or after fetching the pinned dataset and building:
(cd tests/kalman_tfg_shared_modes && ./kalman_tfg_shared_modes-sim --regression)
python3 tools/tfg_shared_modes_sim_table.py
(cd plots/kalman_tfg_shared_modes && bash draw_plots.sh)
cp plots/kalman_tfg_shared_modes/*.pgf doc/kalman_tfg_shared_modes/
```

The pinned v1.3.0 5 kn underway vessel-RAO inputs contain eight scored JONSWAP/PM-Stokes
records, 1200 seconds at 200 Hz. Scores and gates use the final 900 seconds.
The TFG limits in the simulator are its own measured worst errors plus 0.5%,
rounded upward to 0.0001; provenance, per-record scores and limits are in
`reports/results/tfg_shared_modes/deterministic_regression.json`. They cover
vertical and 3-D displacement, roll/pitch/yaw, accelerometer Z and 3-D bias,
and the gyroscope's own 3-D bias error. They are regression sentinels on fixed
records, not statistical population bounds or independent validation results.

`run_tests.sh` verifies the input hashes, requires eight actual passing gates,
checks fail-closed handling of missing/short inputs and altered protocols, and
corrupts one record's displacement truth to demonstrate that numeric gates fail.
No-record and nonfinite-error runs return nonzero. `--regression` rejects sensor,
seed, excitation, variant, lever-arm, tuning and engine overrides that invalidate
the fitted protocol. Seed changes, `--nomag`, `--no-noise`, `W3D_CSM_VARIANT=CSM_OFF`
or `CSM_NOCONF` without `--regression` are explicitly labeled diagnostic runs.

The standalone method article is
`doc/kalman_tfg_shared_modes/kalman_tfg_shared_modes-method.tex`. The original
comparative study remains `kalman_tfg_shared_modes.tex`; the replay gallery is
`kalman_tfg_shared_modes-charts.tex`. CI runs this directory's regression before
regenerating the replay table, figures and all three PDFs.
