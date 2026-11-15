# Moving-vessel replays

All current external-record runs use oceanography-waves-lib **v1.4.0** at 5 kn,
with the 28 ft vessel RAO at the centre of gravity. The two independently
checksummed archives contain the same 30 CSV names but different trajectories:

| Case | Archive | Mean heading | Steering |
| --- | --- | --- | --- |
| `following` | `sim-data-files-vessel-rao-28ft-5kn.zip` | north, 0° | ±3°, 20 s sine |
| `head` | `sim-data-files-vessel-rao-28ft-5kn-south.zip` | south, 180° | ±3°, 20 s sine |

Regular waves in the southbound case are head seas; directional spectra are
centred on bow-quartering seas. The archive's incident height and period are
not the vessel's encounter height and period.

`make all` builds and tests the project, then runs every central replay profile
on **both cases**. `make sim-regression` runs only that matrix. It writes
`reports/results/sim_v1.4.0/{following,head}/`: input and executable provenance,
full-rate metrics in JSON/CSV, seven metric charts, and per-profile logs.
Any missing record, nonfinite metric or failed numeric gate fails the command;
complete failed measurements remain available. Normal CI never refits limits.

For one case, run:

```sh
python3 tools/sim_matrix.py --cases head --build --jobs 2
OCEAN_IMU_SIM_CASE=head python3 tools/simulation_studies.py carma_shared_modes_validation
```

`OCEAN_IMU_SIM_CASE` selects one process's input hashes, reference trajectory
and per-record gates. Set it **before starting Python or a simulator**. The
standalone default is `following`. `tools/sim_dataset.py --case all --dest DIR`
creates `DIR/following` and `DIR/head`; it never extracts both into one folder.
The study and build workflows use separate jobs and case-named artifacts.
The existing following-sea PDF/SVG release filenames stay stable; head-sea
products carry a `-head` suffix.

The scoring reference subtracts only the analytically declared cruise
**displacement and velocity**. All total acceleration and recorded IMU samples
remain available to the estimator. The Python surrogate generator and C++
scoring reference share the 20 s period and the 180° rotation. Neither guesses
the case from identically named CSV files.

The study inventory is in `tools/simulation_studies.py`: paired OU validation
and robustness, TFG comparison, model mismatch, engine noise, startup, wave
direction, low–high–low transitions, lever arms, CARMA shape/excitation/coupling
ablations, shared-mode degradation, paired QMEKF/TFG dataset comparisons, and CSV timing/live-bias audits.
Each is run against both cases by `simulation-studies.yml`. Synthetic-only
prior-process, Lie-group, calibration and stability experiments retain their
independently defined scenarios; they do not consume either ZIP.

Earlier committed study bundles keep their recorded input hashes and remain
historical evidence until replaced by a complete run. They must not be relabelled
as v1.4.0 results. Current case-specific evidence lives under `sim_v1.4.0`.

Historical counterfactual and tuning-search tools retain their explicitly selected controls and configurations. They are not current acceptance studies and are not automatically rerun with substituted inputs.
