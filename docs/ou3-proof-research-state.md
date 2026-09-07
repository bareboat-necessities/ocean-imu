# OU-III proof research state

## Current hypothesis

Complete BRMM replaces spectral SEA3 membership underneath the primary
bounded-bias practical-motion theorem. The full shipping active 21-state
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
Existing SEA3 PASS flags retain their original scope; all BRMM gates are false.
Quiet, oscillatory, mixed windows, both modes and the closed estimate ball
need explicit coverage. Scalar motion recurrence does not imply vector PE.

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
Actual bias/BRMM/retained-SEA3 LaTeX inputs compile in a 14-page syntax smoke
build; this is not a complete manuscript build. GitHub eight-case results
are pending.

## Failure analysis and limiting quantity

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
- Legacy SEA3 linear/Stokes source tools and matrix gates remain scoped
  specializations, not BRMM admission or BRMM P3/P4 proofs.

## Alternatives and next falsifiable experiment

The strongest criticism of the preceding method is that no amount of outward
tightening restores correlations removed from its port domain. Distinct
alternatives are: (1) source-primitive/recurrence equalities and dense joint
sectors; (2) independently weighted channel gains after preserving the
physical graph; (3) a physically structured/path-dependent motion storage.
BRMM supplies a new source premise for (1); it does not make the old
independent-port witness useful.

First obtain the eight-case CI statistics without changing the recordings
or frozen caps. Record actual sampled-envelope violations and unresolved
windows. Establish continuous no-DC/recurrence membership from the generator
or qualified physical source, with same-history primitive/bias/frontend
attachment. Then test the connected motion master with these joint
constraints and explicit channel budgets. Require a useful composed bound
and prefix retention margin before source-uniform enclosure. P4/P5 remain
open; do not promote them from a successful sampled audit.
