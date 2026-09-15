# OU-III ALT current handover

Read `AGENTS.md`, this file, `ou3-alt-deployment-prerequisite.md`,
`ou3-alt-attitude-atlas.md`, `ou3-alt-startup-disturbance-contract.md`,
`ou3-alt-wpe-uniform-supplies.md`, `ou3-alt-startup-pre-rho.md`,
`ou3-alt-proof-plan.md`, and the ALT section of
`ou3-proof-research-state.md`. The original P2/P3/P4/P5 route continues
independently through `ou3-brmm-main-handover.md`.

Continue from a fresh branch off `main`. The research question is completion
of the finite pre-rho qualifications, not estimation of rho. The numerical
startup profile has already been selected and authorized; do not ask the next
conversation to choose it again. This document is the current starting point;
PR #528 contains the implementation and validation record.

## Controlling result

Signed magnetic-counter safety is closed for every finite event prefix. The
shipping attempt, accepted-sample and rejected-sample counts saturate at INT_MAX;
measurement, statistics and release checks continue. The exact induction and
source audit are in `ou3-alt-deployment-prerequisite.md`. No maximum call rate
was added. The master has no falsified counter prerequisite, but universal
startup and source-uniform deployment arithmetic still block rho estimation.

Startup uses the commissioned sensor profiles in
`ou3-alt-startup-disturbance-contract.md`: residual vector caps 0.30 m/s² and
0.02 rad/s for BMI270, 0.50 m/s² and 0.03 rad/s for MPU6886, plus a separate
same-history total direction-error mean/primitive budget 0.10 / 1.5 s.
These are declared engineering requirements supported by datasheet evidence;
hardware admission is not inferred. Every admitted sample has a strictly
positive seed-norm margin. The arbitrary-bounded startup obstruction remains
a regression against conflating startup with post-Live ISS. The existing
Mahony certificate's smaller seed region does not cover these profiles, so
actual capture needs a new proof. The master has 32 closed components and
11 open top-level qualifications.

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

The admitted source factory also accepts a certified empty pre-Live magnetic
prefix. Its shared magnetic event composer carries waiting calls, later north
and saturated counts into the following source-owned IMU event. These are exact
conditional program relations; they do not prove universal startup/capture or
source-uniform arithmetic.

## Retained work

The reset-rooted guard/private-Mahony/LPF/stillness/band/WPE/TuneState product
crosses goLive without new frontend snapshots. Regional Mahony, physical source
restriction, BIAS0/1/2, physical prediction, inverse-free measurements and
rank-three covariance relations remain available. Regional Mahony retains its
original sensor/seed premises; it is not a certificate for the larger
commissioned startup profiles.

The RN32 clock reaches the default timeout comparison at sample 30,002; shared
startup-plus-600-IMU budgets end at 30,602. That deadline is conditional on
source-produced gravity alignment. The magnetic 8.5-uT / 0.63-rad accuracy
calculation remains conditional on a small accumulation-to-handoff frame bound.
That accuracy premise is not needed by the finite atlas word: proper rotations
give the universal chord bound 2 and the same-state mean perturbation bound
157 uT. The graph retains the actual frame and magnetic discrepancy, without
claiming this coarse bound proves nonvanishing north or useful contraction.

## WPE arithmetic qualification

The final tuning-frequency discrepancy is uniformly bounded by the outer
clamps, including nonfinite fallback and retained statistics. For exact bounds
[L,U] and compiled bounds [RN32(L),RN32(U)], the same-history residual lies in
[RN32(L)-U, RN32(U)-L]. This proof requires neither branch agreement nor a libm
accuracy estimate. The two-clamp accepted-update relation additionally uses the
statistics interval. Neither bound proves execution totality or a useful
contraction margin. The master has 11 remaining open qualifications.

The moment graph retains shipping's literal `(alpha*v)*v` and `(alpha*eta)*eta`
order. Each machine history carries its own post-log-update usable latch,
including both inclusive elapsed/history comparisons and the one-way retention
rule. Startup and Live compose independent machine/exact production, initialization
and frequency choices. The strong wrapper supplies its carried machine WPE
state to the lower tuner, then checks the same source-driven moment/log
successor. The eager getter and invalid-result prior fallback are retained.
See `ou3-alt-wpe-branch-composition.md`. Agreement of comparisons is no longer
a prerequisite; source-uniform arithmetic and target qualification remain open.

`ou3-alt-wpe-uniform-supplies.md` supplies a reset-rooted induction for every
bounded vertical input with |x|<=32, including moments, raw period and log/exp
domains. The same sensor packet and literal Mahony normalization/projection
give that input bound after each defined scalar observer update, without a
tilt-accuracy assumption. Preceding observer totality, hardware arithmetic and
the complete source word remain open; the large finite bounds are not ISS gains.

## Exact open qualifications

`tools/stability/ou3_alt_contraction/finite_master_guard.py` is authoritative.
Its current eleven `open_obligations` are:

1. `source_uniform_timeout_aligned_branch_reachability`
2. `near_antiparallel_Eigen_JacobiSVD_solver_correspondence`
3. `universal_startup_source_and_branch_reachability`
4. `startup_source_uniform_deployment_supplies`
5. `target_libm_and_compiler_profile_correspondence`
6. `WPE_machine_target_libm_and_compiler_selection`
7. `WPE_source_uniform_machine_supply_bounds`
8. `WPE_log_and_exp_libm_correspondence`
9. `Qaxis_exp_libm_correspondence`
10. `all_event_arithmetic_witnesses_source_uniform`
11. `complete_source_uniform_600_step_word`

The 32 closed entries are subordinate components. Their count does not imply
that any of these eleven is closed. In particular, conditional WPE induction
does not prove preceding Mahony totality or continuation of the startup sensor
bounds on an arbitrary post-Live ISS word.

## Next falsifiable work

Start with actual startup alignment/capture under the selected sensor and
temporal contract. Before seed arithmetic, the two seed cones are approximately
71.99 and 76.26 degrees. Even retaining the old temporal budgets, the old
Mahony metric needs seed levels 2.05463 / 2.27142, exceeding both its 1.7689
outer level and 1.87267 chart ceiling. Raising that quadratic level alone is
a rejected route; reducing the chosen error caps to fit it is not justified.

Compare a nonlinear observer storage, a same-history finite-time source/LPF
argument and an informative-interval recovery argument against the literal
gravity-alignment predicate. Any recovery premise beyond the selected contract
must be stated and justified, not silently assumed. Follow the failure-analysis
and critic protocol in `AGENTS.md` before refining a failed method.

In the arithmetic branch, attach the actual target/compiler/Eigen/libm profile
and qualify the seed/SVD/nonfinite branches, guard/Racc/frontend/tuner supplies
and complete magnetic history. Preserve the existing scalar and FMA branches,
source/bias ancestry, atlas transport, saturated counters and frontend memory.
Then compose the complete same-history 600-transition word and require the
finite master to pass `assert_finite_storage_master` before storage or rho work.
In charts 1..3 coordinate zero is a 180-degree error, not zero attitude error.
Every-prefix retention, useful ultimate bounds and indefinite no-restart tiling
remain separate later requirements, including wrapper-clock lifetime treatment.

## Reproduction and validation boundary

From the repository root:

```sh
export PYTHONPATH="$PWD:$PWD/tools/stability:$PWD/tests/ou3_alt_contraction"
python3 -m tools.stability.ou3_alt_contraction.finite_master_guard \
  --output /tmp/ou3-alt-pre-rho.json
python3 -m unittest test_finite_startup_sensor_contract \
  test_finite_wpe_uniform_bounds -v
```

The master command should validate its report while keeping all ALT PASS and
storage-readiness flags false. Report validation is not theorem completion.
The exact finite-algebra CI selection is in `.github/workflows/ou3-alt-contraction.yml`.
For native correspondence set `OU3_ALT_REQUIRE_NATIVE=1` and, when necessary,
`EIGEN_INCLUDE_DIR` to the installed Eigen directory. `make all` is the required
repository build; use `CPPFLAGS+=' -I<eigen-directory>'` if Eigen is not on the
default include path. `tools/sim_dataset.py` fetches/verifies the simulation data.

The inherited ALT suite still has phase-1 storage, covariance-envelope,
endpoint-partition and magnetic AD failures. The original live-entry audit
still questions H18-to-A21 release and session-origin S. Keep these gates
intact; finite component tests and host/MCU compilation do not establish target
arithmetic or any of the universal proofs. Exact run results and links belong
in PR #528, rather than being copied into a new historical handover file.

Zero wind heel, dormant guard, zero lever arm, full-21 covariance, joint24
motion/bias information, actual-applied R_S, corrected COMPLETE-BRMM and frozen
P3 delta=1e-18 remain unchanged. All ALT PASS and storage-readiness gates are
false. Validation results belong to the active PR description and CI logs.
