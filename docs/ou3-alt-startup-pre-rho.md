# ALT startup obligations before rho

The single-Cayley entry obstruction is resolved by the four-chart finite
runtime in `ou3-alt-attitude-atlas.md`, including ungauged Live and later north
acquisition. The obstruction below explains why one chart is insufficient.
`assert_finite_storage_master` remains blocked on universal source/deployment
qualification. The original P2/P3/P4/P5 track is independent and unchanged.

## An admitted physical family that reaches the chart pole

Take a level boat at rest with any constant true heading. For all physical time,
let `p=v=a=0`, gyro bias and accelerometer bias be zero, and sensor disturbances
be zero before conversion to the declared machine format. The wave displacement
and every centered primitive are zero. This satisfies the quiet-zero member of
corrected COMPLETE-BRMM and the zero solution of BIAS0, BIAS1 and BIAS2; there is
no nonzero DC displacement or reset of S. Rate is zero, wind heel and lever arm
are zero, and the constant accelerometer stream leaves the vibration guard
dormant. A fixed world magnetic field `(15,0,20) uT` satisfies MAG-BMM150-DET-v1.

MAG-CALL-SCHEDULE-v1 imposes post-gauged-Live service, not a pre-Live acquisition
deadline. A startup history with no magnetic calls before timeout is therefore permitted
by the current source language. Calls can subsequently continue at ordinary 25 Hz;
the obstruction requires neither a huge call count nor an artificial burst.

Every heading gives the same pre-Live sensor history: gyro zero and specific
force `(0,0,-g)`. Ordinary FromTwoVectors seeds the private observer with identity.
Its correction cross product and integral remain zero; its quaternion remains
a nonzero scalar quaternion. In the named scalar binary32 graph, after the
second update it is exactly

`q_proxy=(16748919/16777216,0,0,0)`.

The next update reproduces that q, integral and vertical output exactly. Those
recurrences do not read elapsed time, so this fixed-point identity proves their
continuation without enumerating a replay. The normalized public proxy is
identity. The gravity-aligned branch holds, while no yaw gauge has been learned.
The timeout permits handoff independently of tuner readiness and north.

The unchanged public wrapper regression executes construction through this
handoff, without installing a state. Under scalar Eigen 3.4 and no FP contraction
it reaches Live at sample 30,002 with identity quaternion and no north lock.
This host correspondence supports the program relation; target ESP32/Eigen/libm
qualification remains a separate obligation.

Choose true world-to-body quaternion `(0,0,0,1)`, a perfectly valid 180-degree
yaw with zero tilt. The fresh error quaternion then has scalar component zero.
The local `CORE.cayley` operation defines

`c = 2 q_error.vector / q_error.scalar`.

Its denominator is exactly zero. This is not a degenerate physical quaternion.
Moreover the unit rational family

`q_n=(2n/(n²+1),0,0,(n²-1)/(n²+1)), n>=1`

has the same constant physical IMU history and `c_z=n-1/n`. Thus deleting only
the exact pole would still leave no uniform finite fresh Cayley radius. A
smaller assumed entry set would change the theorem source.

## Literal timeout and finite clocks

Shipping uses `max(proxy_timeout, settle + 2*max(mag_min_window,1) + fallback)`
when magnetometer service is enabled. Defaults give an acquisition floor of
60 seconds and total threshold 150 seconds. Timeout additionally requires an
initialized proxy and the source-produced gravity-aligned branch. It requires
neither TunerReady nor north. The conditional bridges carry a control decision
for Cold/TunerWarm handoff and retain the same frontend, pending TuneState,
WPE/band/statistics and a_w-sync ancestry.

For `t[k+1]=RN32(t[k]+RN32(0.005))`, exact enumeration gives

| Quantity | Exact value |
| --- | --- |
| Last clock below 150 | `t[30001]=9830121/65536` |
| First clock at least 150 | `t[30002]=9830449/65536` |
| Physical time at first crossing | `150.01 s` |
| End of the following 600 IMU transitions | sample `30602`, physical `153.01 s` |

All finite machine-history budgets use this common conditional horizon. The
tau induction remains strictly inside its positive local domain after the two
additional updates. None of this proves the aligned predicate for every source
at the deadline. Later handoffs need either a separately proved deadline or
longer/indefinite arithmetic supplies.

## Magnetic accuracy needs the full frame history

For an accepted magnetic sample let
`A_i=R_tilt_hat_i R_true_WB_i` and let `G_L=Rz(-psi_true_at_handoff)`.
The valid conditional bound is

`||A_i-G_L||_2 <= 2s  =>  ||mean-G_L B|| <= 2 B_max s + b_HI_max + n_max`.

It retains heading variation during accumulation and the delay to handoff.
The declared 0.02-rad gravity-direction error does not supply this full-frame
premise. With perfect tilt, north and south samples of B=(15,0,20) have mean
(0,0,20): tilt alone cannot ensure a nonvanishing mean horizontal field. This
two-rotation algebra refutes the premise substitution; it is not presented as
a complete shipping startup trajectory.

The arithmetic implication from `s<=0.01` and handoff tilt `<=0.02` still gives
8.5 uT perturbation, yaw `<0.61 rad` and full angle `<0.63 rad < pi/4`. The API
requires those full-frame/handoff bounds explicitly and no longer labels them
source-qualified consequences of the gravity bound.

The finite atlas word does not need this small-angle accuracy premise.
For every pair of proper rotations, ||A_i-G_L||_2<=2. The same arithmetic mean
therefore has the unrestricted perturbation bound 157 uT. This bound need not
keep north nonzero, but both the not-ready branch and every nonzero fresh
attitude are represented. The actual frame and reference discrepancy remain
in the graph. Small-angle magnetic accuracy may matter for a later basin or
usefulness theorem; it is no longer an independent pre-rho master prerequisite.

## Remaining work, in dependency order

1. The four-chart cover and conditional ungauged event continuation are
   attached. Keep their chart index and exact transport in every remaining
   source/machine layer; do not replace chart-dependent physical attitude by
   coordinate magnitude alone.
2. Use the declared commissioned startup profiles and proved seed-norm margins
   in `ou3-alt-startup-disturbance-contract.md`. Their larger seed set is not
   covered by the original conditional Mahony invariant. Prove actual
   alignment/predicate production, seed/SVD/nonfinite branches and finite
   reachability without tightening the chosen sensor caps to fit that invariant.
3. Qualify all target/compiler arithmetic and supplies: WPE exp/log/sqrt,
   band/statistics/tuner, Q-axis, Racc/guard displacement, clocks/scheduler,
   trigonometry/normalization, Eigen LDLT/eigensolvers/floors, and comparisons.
   The bounded-input WPE induction and conditional source projection are in
   `ou3-alt-wpe-uniform-supplies.md`; they leave upstream totality and target
   correspondence open. There are 11 top-level prerequisites remaining.
4. Counter safety is closed by the source-bound saturation induction in
   `ou3-alt-deployment-prerequisite.md`, including accepted/rejected acquisition
   and refinement counts. Retain its projection from mathematical event count
   and compose every literal event on one physical, bias, disturbance,
   covariance and frontend history. The certified empty startup prefix and
   subsequent north acquisition are attached to the admitted source word.
5. Only after the complete master passes the guard may common joint24 storage
   and rho feasibility begin. Every-prefix retention, useful ultimate bounds,
   and indefinite no-restart tiling follow separately.

The pole and nearby family remain regressions against a universal single-Cayley
fresh-entry bound. The current atlas covers both. No ALT PASS is claimed.

rho has since been measured ahead of this dependency list, and the declared
joint24 contraction is falsified on the admitted ungauged word; see
`ou3-alt-word-rho-feasibility.md`. The obligations above stay correctly stated
but are subordinate to `rho_w`, so completing them is no longer the controlling
work.
