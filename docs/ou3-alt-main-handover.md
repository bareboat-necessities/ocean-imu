# OU-III ALT current handover

Read `AGENTS.md`, this file, `ou3-alt-attitude-atlas.md`,
`ou3-alt-startup-pre-rho.md`, `ou3-alt-proof-plan.md`, and the ALT section of
`ou3-proof-research-state.md`. The original P2/P3/P4/P5 route continues
independently through `ou3-brmm-main-handover.md`.

## Controlling result

The finite runtime now represents every nonzero attitude error in a four-chart
joint24 atlas. Maximum-component selection bounds each attitude coordinate by
2 and covers the admitted south-heading timeout entry. Prediction, accepted
measurement, tilt reset, magnetic yaw rewrite and H18/A21/hold transitions carry
the chart correctly. All motion/bias information and full 21-state shipping
covariance are retained. The local Cayley obstruction remains a regression,
not an unresolved global representation problem.

The interleaver now carries ungauged Live and delayed initial north acquisition.
Initial acquisition and the gravity gate after Live use MEKF tilt; continuous
calibration and refinement use private-observer tilt. The same packet can
establish north, refine, apply calibration and reach the inner MEKF without
accumulating continuous statistics twice. The north-service clock starts at the
actual gauge event without moving the physical Live/S origin.

These are exact conditional program relations. They do not prove universal
startup/capture, source-uniform arithmetic or counter safety. Rho remains blocked.

## Retained work

The reset-rooted guard/private-Mahony/LPF/stillness/band/WPE/TuneState product
crosses goLive without new frontend snapshots. Regional Mahony, physical source
restriction, BIAS0/1/2, physical prediction, inverse-free measurements and
rank-three covariance relations remain available.

The RN32 clock reaches the default timeout comparison at sample 30,002; shared
startup-plus-600-IMU budgets end at 30,602. That deadline is conditional on
source-produced gravity alignment. The magnetic 8.5-uT / 0.63-rad accuracy
calculation still requires a full accumulation-to-handoff frame bound.

The signed magnetic counter graph rejects INT_MAX increment instead of silently
extending C++ with Python integers. Avoidance of that edge is unproved under
MAG-CALL-SCHEDULE-v1, which has no uniform upper count.

## Continuation

1. Prove startup source/control reachability, including seed/SVD/nonfinite
   branches, and qualify complete initial/refined/continuous magnetic history.
2. Attach the actual target/compiler/Eigen/libm and caller profile; close
   source-uniform arithmetic supplies, comparison branches and finite counter
   safety. A sketch's one-mag-call-per-IMU bound cannot silently replace the
   current broader asynchronous language.
3. Compose the entire same-history 600-transition machine word and require
   `finite_master_guard.build()` / `assert_finite_storage_master` to pass.
4. Only then construct coercive storage respecting chart semantics and exact
   chart transport. In charts 1..3 coordinate zero is NOT zero attitude error;
   treating every chart origin as zero in one quadratic would be invalid.
5. Every-prefix retention, ultimate bounds and indefinite no-restart tiling
   follow; they still require clock/counter lifetime treatment.

Zero wind heel, dormant guard, zero lever arm, full-21 covariance, joint24
motion/bias information, actual-applied R_S, corrected COMPLETE-BRMM and frozen
P3 delta=1e-18 remain unchanged. All ALT PASS and storage-readiness gates are
false. Validation results belong to the active PR description and CI logs.
