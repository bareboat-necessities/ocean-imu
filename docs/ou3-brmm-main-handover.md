# BRMM stability proof handover

## Starting point

The conditional certificate below was rebuilt on main
85e040371812d1fa32508d521acded0a6328102f. Read this handover and the research
ledger, then check the exact revision and its CI before continuing.
The runtime audit qualifies coverage of the reference simulations; it does
not establish all Normal-Live premises or complete the nonlinear theorem.

## The theorem being pursued

BRMM is the only primary physical motion hypothesis. No retired source
identifier remains in tracked text, filenames or workflows; a regression
test guards this. Spectral/linear-vessel/Stokes modules are reference models,
not mandatory BRMM membership conditions.

The target is bounded accelerometer-bias error plus practical ISS of the
other 18 motion/attitude/gyro-bias errors. The active estimator remains full
21 state. Keep all covariance/gain cross terms, every due S correction with
actual anisotropic R_S, physical -S_true forcing, finite attitude resets,
full process covariance/floors, and the closed 0.4 estimate projection ball.
No bias-error decay and no zero model-mismatch/noise floor are required.
The observed 2--3 percent residual is not a certified uniform accuracy bound.

BRMM-0 requires a_m=dv_m/dt and uniformly bounded same-history v_m over the
entire continuation. It prohibits fixed physical acceleration DC in Q, O
and mixed histories. Constant measurement offset remains BIAS0. Small AC
energy and finite-window impulse bounds alone do not prohibit DC. Do not
subtract a sample mean to manufacture membership. Bounded velocity alone
does not bound position or S; their shared primitives/budgets remain needed.

The current conditional P3 bounds are acceleration 4 m/s^2, body rate 30 deg/s, positive
H_s<=8.5 m plus calm motion. Ten seconds, RMS .03/.05 and impulse 2 m/s
are audit candidates; recurrence/lobe/primitive constants remain unfrozen.
The shipping tuner .03--1.2 Hz is not a physical motion bandlimit.
BIAS0/1 qualification measurements gate deployment, not conditional math.
BIAS2 must be proved on the actual nonlinear history if used as a sector.

## P3's exact scope

Read `ou3_brmm_p3_premises.py`, the primary BRMM LaTeX section and
`ou3_brmm_riccati_metric_p3.py`. The implication is:

BRMM + explicit Normal-Live execution premises implies
Omega_W - 1e-18 P_W >= 0 for the full H18/A21 covariance word.

Normal-Live is a theorem regime, not merely a runtime Live flag. Premises
include accepted accelerometer samples, asynchronous accepted-vector PE with
the declared geometry and recurrence, transparent vibration guard/zero lever
arm, shipping committed parameter/scheduler invariants, successful updates,
all actual R_S/Q/floor/reset events, source-generated entry and exact release.
BRMM scalar recurrence does not imply vector PE. Geometry bounds must attach
to actual measurement Jacobians on the same history. P3's conclusion and
bias-error decay are not premises.

The quantitative route uses the full H18 prior-free 18x18 outward LDLT after
the complete 3-second information word and following prediction; the shipping
hold then permits the exact H-to-A 21x21 covariance bridge, and subsequent
event algebra preserves the margin. Delta remains 1e-18. No finite source
enumeration or spectral response membership is needed for this implication.
Q/O/mixed coverage is conditional coverage, not physical admission.

The exact shipping bias projection changes nominal bias and leaves covariance
unchanged. The formal covariance-word comparison survives the closed boundary;
its transport matrix is not the derivative of the projected nonlinear error
map. P4 still owes that map's storage attachment and all-prefix retention.
Do not present this covariance fact as nonlinear error contraction.

## Existing evidence and what it does not show

The eight v1.1.3 reference records contain 240000 samples each and 1904008
sample-aligned 10-second windows total. All windows are O at both thresholds;
all have a positive explicit opposing-lobe witness. There is no Q coverage.
Minimum sampled outward energy lower bound is 1.4596516559 m^2/s^3 and minimum
J lower bound is .2472575910 m/s at diagnostic chi=1. These are CSV-point ZOH
bounds, not physical intersample or infinite-time source certificates.

Six whole records exceed acceleration 4, all eight exceed body rate 30, and
six exceed candidate 10-second impulse 2. Largest sampled values round upward
to 16.2162 m/s^2, 198.355 deg/s and 8.0912 m/s. This raw source audit does
not attach runtime state. The separate runtime audit below now shows that
violations also occur in Live. These remain admission failures, not
admitted-word counterexamples or filter instability. Keep all violations.
Per-case evidence and CSV hashes are in `docs/ou3-brmm-reference-audit.json`.
Audit run 34169896915 artifact 10035405278 has SHA256
`e4a95423a06e16e6cead32b3e73bd3cee6a579dd2c72d89a3717c6cbca3b5b18`.

The preceding connected point-coefficient supply test from #502 is feasible
but quantitatively unusable: factors .9999569976489486/.9788191291615017 and
common gains 4294967296/68719476736 imply candidate storage bounds at least
1.4037e15/2.2313e13 on the observed forcing. These are not actual filter floors
or certified nonlinear contraction factors. Independent input-energy ports
discarded coupling between latent increments and S_true; one common gain
charged every channel at the most expensive scale. No interval refinement or
metric grid is justified around that witness. Actual R_S corrections were
retained; the previously missing forcing attachment was -S_true.

## Current continuation

Read docs/ou3-proof-research-state.md and docs/ou3-brmm-runtime-audit.json.
The eight-case runtime attachment is now implemented. Live itself contains
acceleration/rate violations, and actual replay Racc/Rmag differ from the
configured P3 values. Do not repeat the raw-window lobe search or select
samples using the very caps being tested. The next question is physically
qualified, separate before-Live/Live envelopes and the corresponding actual
measurement configuration. The large JONSWAP peak includes a second-order
correction larger than the linear contribution and needs model qualification.
A 250 deg/s candidate fails the existing two-occurrence PE transport bound;
full accepted-vector transport is needed before claiming a broader P3.

## Verified P3 and outstanding validation

Exact main 85e040371812d1fa32508d521acded0a6328102f passed conditional BRMM
P3 in run 34173418232 / job 101899163631, including 56 tests. H18 and A21
delta are both 1e-18. H18 worst outward LDLT pivot is 4.987499868870966e-14;
first active A21 bias margin is 1.2499987189052501e-9. Artifact 10036647100
was downloaded and its ZIP SHA256 verified as
03d4011e49f0c54224cc129709631a26d5ee8b5399a435e9e083f22026b8ab03.
This supersedes the queued P3 status at landing. The result is conditional
on its declared configuration and is not a broader-domain or nonlinear PASS.

The separate inherited replay-provenance failure requires genuine full
validation/robustness regeneration through the existing full-study workflow.
No quality gate, dependency hash or proof threshold may be bypassed. See
the current research ledger for the runtime experiment and its validation.
