# BRMM stability proof handover

## Landing and starting point

Work is in [PR #504](https://github.com/bareboat-necessities/ocean-imu/pull/504),
stacked on [PR #502](https://github.com/bareboat-necessities/ocean-imu/pull/502).
Both are separate from merged #500. Main at the last sync is
`f5dacecb7fe73e6daf1e3caeb1b9124d5511c966`; #502 head is
`998eb8c2737f159728a6804e43419a3254b9f8dd`. Neither PR has been merged by this
work. Merge preparation is not a claim that main already contains the changes.

After checks and merge authorization, land #502 first, retarget #504 to main,
recheck its diff and CI, then land #504. A new thread should fetch main and
verify it contains this document and `tools/stability/ou3_brmm_p3_premises.py`
before continuing. If either is absent, the stack has not fully landed.
Do not start from an older theorem or silently discard the stack's changes.

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

Physical bounds are acceleration 4 m/s^2, body rate 30 deg/s, positive
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
to 16.2162 m/s^2, 198.355 deg/s and 8.0912 m/s. No runtime Normal-Live predicate
was attached, so these are source-envelope failures over whole raw records,
not admitted-word counterexamples or filter instability. Keep all violations.
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

## Next decisive work

1. Attach the actual runtime Normal-Live predicate to the unchanged eight
   histories. Record outer/inner Live, held/active bias, actual guard action,
   accepted vector events, hard re-lock/reset events and every S correction.
   Expose cap-exceeding intervals; do not prune them to claim admission.
2. Establish continuous no-DC/recurrence membership and common physical
   primitives. A fixed harmonic model's all-time velocity bound B1 helps,
   but uniform root bounds and exact physical/front-end attachment are needed.
3. Feed those primitive equalities, BIAS0/1 history and any proved BIAS2 sectors
   into the full connected nonlinear motion master. Retain -S_true and actual
   P/H/R/K/reset/projection from every selector on the same lineage.
4. Test useful endpoint motion factor/channel gains and a resulting practical
   floor compatible with chart/domain retention before outward source covering.
   Every prefix needs finite gain and retention, not contraction. P4/P5 stay
   closed until those obligations are certified uniformly over admitted leaves.

Follow AGENTS.md: record exact failures and what they invalidate, perform a
critic pass, and respect the two-strike rule. Do not retune the filter, weaken
quality gates or infer sensor qualification from configured tau_b=5000 s.

## CI and merge gates

The migration revision is `0db18980f4ff7e018e392a847ed4ba0d6ad24ef7`, tree
`de271057b08b870de9fb05f4099fc8518d7f3df3`; GitHub run 34172553221 rebuilds
the theorem/source/P3 chain. Results and any following fixes are recorded in
the research ledger and PR body; do not transfer a PASS from another tree.

Mandatory local `make all` fails because Eigen/Dense is unavailable. CI owns
the C++ build result. The initial migration's Python quality job identified
one unused reference-only local; remove it without changing the mathematics.

An independent OU validation gate reports stale replay provenance for OU-II,
OU-III, the tuner and the OU-III simulation driver. It requires a genuine full
validation/robustness regeneration. The existing `ou-full-evidence-branch.yml`
workflow is opt-in and runs about 1150 simulator replays; main has its automatic
full-study path. Never repair this gate by restamping hashes or bypassing the
evidence contract. Until the required checks are green, the stack is not
merge-ready even if conditional P3 passes.
