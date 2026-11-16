# Steering-separated replay limits

The limits are locked to the completed native replays of the steering-separated runtime,
using the existing 0.5% allowance rounded upward to 1e-5. All 544 cells are retained;
272 horizontal cells additionally gate geographic wave Y RMS directly in metres.
Every one of the 272 per-record 3-D displacement caps is tighter.

This is an explicit baseline change. It does not fix the remaining slow bias or
attitude errors. In particular, 213 accelerometer 3-D bias caps and 200 pitch caps
increase. Their largest relative increases are 60.4% and 62.8%, respectively.
The before/after estimator measurements remain in `../heading_kinematics/runtime.json`.
`heading_limit_update.json` records the cap-change counts and largest regressions.
The qualifying `baseline_*.json` receipts retain their original failed-gate outcomes;
fresh enforced runs, with actual exit codes, are published in `../sim_v1.4.0/`.

## Global default-profile 3-D displacement ceilings

Values are percent of the maximum wave-reference 3-D displacement. Each ceiling
covers both directions; the per-record cap is usually much lower. Interval profiles
are also regenerated in `CurrentReplayLimits.h` and `limits.json`.

| Family | JONSWAP before | JONSWAP after | PM/Stokes before | PM/Stokes after |
| --- | ---: | ---: | ---: | ---: |
| carma | 32.59217 | 6.66545 | 33.85218 | 7.21956 |
| csm | 32.15921 | 6.26970 | 33.47564 | 6.66302 |
| ou2 | 71.98221 | 17.35218 | 72.30455 | 17.18685 |
| ou3 | 69.43741 | 20.57078 | 71.03592 | 20.37571 |
| tangent | 32.15172 | 6.26987 | 33.46893 | 6.66346 |
| tcarma | 32.59298 | 6.66569 | 33.85270 | 7.21972 |
| tcsm | 32.15176 | 6.26983 | 33.46844 | 6.66334 |
| tfg | 136.82147 | 18.32304 | 129.39057 | 17.97445 |
| tfg_legacy | 47.95832 | 31.03596 | 47.65723 | 31.03184 |

## Independent Y caps for the 0.27 m JONSWAP record

Production defaults, metres RMS, trailing 900 s. The new Y gate is independent of
the 3-D gate, so an improvement in another axis cannot hide a renewed Y error.

| Family | Following | Head |
| --- | ---: | ---: |
| ou2 | 0.01256 | 0.00630 |
| ou3 | 0.01283 | 0.00832 |
| tfg | 0.01689 | 0.00836 |
| carma | 0.00458 | 0.00236 |
| csm | 0.00431 | 0.00231 |
| tcarma | 0.00458 | 0.00236 |
| tcsm | 0.00431 | 0.00231 |
| tangent | 0.00431 | 0.00231 |
| tfg_legacy | 0.03935 | 0.02583 |

## Per-record changes

Counts include both directions and all applicable profiles. Unchanged non-horizontal
methods remain in their existing inventory; small changes in other native targets
are included rather than hidden.

| Metric | Tightened | Increased | Unchanged | New |
| --- | ---: | ---: | ---: | ---: |
| roll | 231 | 42 | 15 | 0 |
| pitch | 72 | 200 | 16 | 0 |
| yaw | 142 | 130 | 32 | 0 |
| heave | 181 | 103 | 244 | 0 |
| displacement | 272 | 0 | 0 | 0 |
| acc_z | 135 | 135 | 2 | 0 |
| acc_3d | 59 | 213 | 0 | 0 |
| gyro_3d | 253 | 18 | 1 | 0 |
| displacement_y | 0 | 0 | 0 | 272 |

## Reproduction

The old tables can be read from commit `e0eeb84de8a96c8b2ec8742069b915646d36f1d2`.
To reproduce the updated tables from the two immutable qualifying receipts:

```sh
mkdir -p /tmp/heading-limit-inputs/following /tmp/heading-limit-inputs/head
cp reports/results/filter_sim_gates/baseline_following.json /tmp/heading-limit-inputs/following/measurements.json
cp reports/results/filter_sim_gates/baseline_head.json /tmp/heading-limit-inputs/head/measurements.json
python3 tools/lock_simulation_baseline.py --measurements /tmp/heading-limit-inputs \
  --write-reviewed-baseline --reason "Reproduce PR #763 steering-separated native baseline"
python3 tools/generate_record_replay_limits.py --check
make all -j1 SIM_JOBS=4
```

Normal CI reads these limits and never invokes the lock command. The source, input
and binary hashes in the receipts distinguish the qualifying measurements from
the subsequent runs that enforce the new gates.
