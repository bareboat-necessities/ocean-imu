# OU-III proof research state

## Current hypothesis

Complete BRMM is the physical source underneath the primary bounded-bias
practical-motion theorem; spectral models are reference specializations. The full shipping active 21-state
filter, covariance/gain cross terms, finite resets/projection, every due
actual anisotropic R_S correction and physical -S_true input remain.
The performance storage measures the other 18 errors. Bias-error decay and
a zero mismatch/noise floor are not requirements.

BRMM-0 requires a_m=dv_m/dt and uniformly bounded same-history velocity on
the whole admitted continuation. It excludes fixed nonzero acceleration DC
in Q, O and their mixtures: every long-time mean is bounded by 2*V_m/T.
Constant measurement offset remains a BIAS0 coordinate. A small AC energy
or fixed 10-second impulse cap alone does not exclude DC. Bounded velocity
alone does not bound displacement/S; retain their shared graph and forcing
budgets. No sample-mean subtraction can manufacture physical membership.

The declared acceleration/body-rate caps remain 4 m/s^2 and 30 deg/s.
Positive H_s<=8.5 m and calm motion are admitted. T_R=10 s, quiet RMS .03/.05
and impulse cap 2 m/s are initial audit candidates; recurrence/lobe/primitive
constants are unfrozen. Shipping tuner defaults remain .03--1.2 Hz, not a
physical motion bandlimit. A .02 Hz sinusoid need not return inside 10 s.

## Evidence and current experiment

The BRMM contract, LaTeX theorem and proof plan separate source admission,
P3 at delta=1e-18, endpoint motion gains, every-prefix gain/retention and P5.
P3 is being rebuilt as BRMM plus explicit execution premises implies the
full matrix inequality, with Q/O/mixed conditional coverage and a checked
covariance identity at the closed bias projection boundary. Source admission,
P4 and P5 remain open. Scalar motion recurrence does not imply vector PE.
The runtime Live flag does not establish the full Normal-Live theorem regime.

The focused BRMM CI audits every sample-aligned 10-second window of all eight
unchanged v1.1.3 JONSWAP/PM-Stokes reference CSVs. It records energy, impulse,
TV, explicit opposing-lobe witness offsets, cap violations and primitive
diagnostics. Outward arithmetic encloses binary64 CSV-point ZOH statistics;
a missing finite-search lobe witness is unresolved, not nonexistence.
This is not physical intersample inclusion, all continuous offsets,
Normal-Live/frontend/tuner admission, infinite-time no-DC, or a source cover.
Empirical margins do not set physical theorem constants.

Eight focused local tests pass, including independent rational-arithmetic
energy/impulse/TV checks, verification of each positive lobe witness, calm
motion, the quiet-DC loophole, and the long-period one-sided window.
Actual bias/BRMM/reference-model LaTeX inputs compile in a 14-page syntax smoke
build; this is not a complete manuscript build.

At c1046c1f, GitHub audit 34169268266 / job 101886317643 succeeds. All eight
records have 240000 samples and 238001 windows each: 1904008 windows total.
Every window is O at both .03/.05 RMS thresholds, with no ambiguous
classification and no missing opposing-lobe witness. Thus these records
provide no Q coverage. The minimum sampled energy lower bound is
1.4596516559 m^2/s^3, and the minimum J witness lower bound is
.2472575910 m/s at diagnostic chi=1. They are ZOH bounds, not continuous
source constants. The full per-case summaries and CSV hashes are retained
in docs/ou3-brmm-reference-audit.json. Every-window NPZ evidence is artifact
10035278812, SHA256
8be76bd010948f36f2ea564fea2dd2c3b0cd5bb814163c2b503e759c18df8409.
The theorem CI job 101886317504 also succeeds (14-page syntax smoke).

| Family | H_s m | Max acceleration m/s^2 | Max body rate deg/s | Max 10 s impulse m/s |
| --- | ---: | ---: | ---: | ---: |
| JONSWAP | .27 | 2.8321 | 36.599 | 1.0083 |
| JONSWAP | 1.5 | 6.7233 | 83.971 | 2.7156 |
| JONSWAP | 4 | 9.2272 | 129.129 | 5.2838 |
| JONSWAP | 8.5 | 16.2162 | 198.355 | 8.0912 |
| PM-Stokes | .27 | 2.1850 | 46.697 | .9870 |
| PM-Stokes | 1.5 | 5.1905 | 123.326 | 2.9807 |
| PM-Stokes | 4 | 6.6333 | 149.547 | 4.8231 |
| PM-Stokes | 8.5 | 8.5683 | 197.461 | 6.8645 |

Displayed maxima round upward; the JSON retains detailed outward sample
bounds and radian-rate bounds. Physical intersample extrema remain unknown.

## BRMM naming and P3 premise migration

The current tree renames the retired physical-source identifiers, imports,
paths and workflow/artifact names to BRMM. Reference spectral artifacts remain
explicit reference-model metadata; they do not define primary BRMM membership.
The canonical source carries bounded motion primitives and one common history.
The new P3 premise manifest distinguishes code identities from execution
assumptions and checks the exact shipping projection body: it changes only
nominal bias, so the covariance comparison is unchanged. The covariance-word
transport is not the derivative of the nonlinear projected error map; P4
must still attach that map and prove retention. The numeric P3 chain is being
rebuilt in GitHub CI, retaining delta=1e-18 and Q/O/mixed conditional coverage.
At 0db18980, source-foundation and the actual theorem syntax build pass.
Quality job 101895520177 reports one F841 unused `sea` local left after
removing reference-only data from P3 mandatory premises. Class: implementation
cleanup defect, not theorem failure. Remove that unused assignment and rebuild;
no numerical premise or proof gate changes. A focused check also caught a
reference hard-set paper-parity marker still expecting the primary-source
notation after the reference section was relabelled. Bind it to the actual
reference hard set and its explicit non-membership scope; retain compactness,
outward representation and rejection of stochastic/seeded surrogates.
Class: documentation-to-contract binding defect, not a numerical failure.

Local verification: 15 focused BRMM/source/premise tests pass, including a
mutation that inserts a covariance write into projection and is rejected.
All proof/test Python files parse. Mandatory `make all` fails at
`src/ahrs/KalmanQMEKF.h:30`, missing `Eigen/Dense` despite the existing
`-I/usr/include/eigen3`. Class: local dependency/infrastructure failure; it
says nothing about theorem feasibility. CI installs Eigen and must provide
build results. No include path, production filter or numerical gate changed.
The existing OU validation replay-provenance mismatch is a separate red gate;
only genuine full replay regeneration may replace that evidence.

## Failure analysis and limiting quantity

The sampled numerical-envelope hypothesis fails: six whole records contain
definite acceleration >4 m/s^2, all eight contain body rate >30 deg/s, and
six have 10 s impulse >2 m/s. For example JONSWAP H_s=8.5 has 30314/48304
acceleration/rate exceedances and 134408 impulse-violating windows, starting
at sample 0 for impulse. Failure class: proposed source-envelope admission,
not an enclosure/conditioning failure. It invalidates the claim that the
eight full reference histories already fit the proposed numerical caps.
It does not falsify no-DC, the lobe condition on these sampled windows,
conditional BRMM stability or the filter. No runtime Normal-Live predicate
was applied, so this is not yet an admitted-word counterexample.

Critic pass: sea height alone cannot guarantee the acceleration/rate caps,
and a 2 m/s impulse cap is quantitatively too small for these records.
Shrinking sample boxes or sharpening rounding cannot resolve the measured
gaps. Keep the frozen caps and do not silently discard offending samples.
Alternatives are (1) attach the actual execution predicate and determine
which source intervals are inside/outside the declared theorem domain;
(2) separately qualify a physically justified broader domain, rebuilding its
dependent gravity/chart/P3/P4 bounds; (3) retain the original domain and
state explicitly which simulation portions it does not cover. The current
next experiment is (1). For recurrence alone, 10.2 m/s is a sampled-data
V_R candidate with 25 percent margin; it is not a frozen physical constant.

The no-DC premise has a separate analytical model route: the retained
fixed-root harmonic graph has a_m=dv_o/dt and ||v_o||<=B1. Constant Stokes
velocity drift contributes no acceleration DC. A uniform root-family B1
bound and exact physical attachment still have to be supplied.

The preceding point-coefficient supply test is mathematically feasible but
fails as a useful accuracy/retention witness. At code head 7898eb06,
connected CI 34165890828 passed 42 tests and both shipping captures.
H18/A21 candidate factors are .9999569976489486/.9788191291615017, with
common normalized gains 4294967296/68719476736. Their composed storage
bounds would be at least 1.4037e15/2.2313e13 using the observed source energy,
even with zero true-bias root. These are unusably loose bounds, not actual
filter error floors or certified nonlinear contraction factors.
Artifact 10034187923 has SHA256
b4bba842068c52f59049cf0cfb4416d2601fc6413293674c2ebd9d38c31c2b28.

Failure class: proof-method/gain-budget failure. Arbitrary independent energy
ports discard physical coupling between latent increments and S_true;
one common gain charges every channel at the most expensive scale.
This invalidates using that witness for a useful uniform funnel. It does not
falsify the filter, bounded-bias theorem, or the still-unproved BRMM target.
No interval refinement or common-gain search is justified for that witness.

BRMM's current limiter is a qualified no-DC/recurrence source graph and its
usable joint forcing/gain budget, with separate Q/P3 coverage. No nonlinear
gain is inferred merely from an impulse statistic or optional BIAS2 sector.

Mandatory local build: make all fails in tests/ahrs while compiling
ahrs-qmekf-sim.cpp, at src/ahrs/KalmanQMEKF.h:30:10:
fatal error: Eigen/Dense: No such file or directory.
Failure class: local infrastructure; preserve Eigen/include paths.
It does not invalidate the math or Python/LaTeX checks. CI installs Eigen.
The separate inherited OU replay-provenance mismatch is not weakened.

## Retained facts and dead ends

- All 137/108 actual R_S corrections are present in the connected H18/A21
  words. A21 net signed S energy is -1.121508422 on that point. The older
  physical attachment omission was the -S_true input, not the R_S gain.
- Actual forced W18 endpoint ratios 7.971765/3.586608 are not homogeneous
  contraction tests. The reset-deleted stronger A21 expander is not a
  canonical counterexample because finite reset/source attachment is missing.
- Qualified true-bias bound plus the .4 estimate clamp gives bounded error.
  The old .35 interior does not cover the primary closed projection ball.
  Assembled-sensor measurements gate deployment, not conditional mathematics.
- BIAS2 may sharpen a gain only through a proved sector on the same corrected
  nonlinear/source history. Positive source AC energy alone is insufficient.
- Legacy BRMM linear/Stokes source tools and matrix gates remain scoped
  specializations, not BRMM admission or BRMM P3/P4 proofs.

## Alternatives and next falsifiable experiment

The strongest criticism of the preceding method is that no amount of outward
tightening restores correlations removed from its port domain. Distinct
alternatives are: (1) source-primitive/recurrence equalities and dense joint
sectors; (2) independently weighted channel gains after preserving the
physical graph; (3) a physically structured/path-dependent motion storage.
BRMM supplies a new source premise for (1); it does not make the old
independent-port witness useful.

The eight-case sampled lobe audit is complete; no longer search for a return
lobe on those same windows. Attach the actual Normal-Live execution predicate
to the unchanged records and expose every cap-exceeding interval, without
selecting an easier source or treating missing admission as a P4 failure.
Establish continuous no-DC/recurrence membership from the generator
or qualified physical source, with same-history primitive/bias/frontend
attachment. Then test the connected motion master with these joint
constraints and explicit channel budgets. Require a useful composed bound
and prefix retention margin before source-uniform enclosure. P4/P5 remain
open; do not promote them from a successful sampled audit.
