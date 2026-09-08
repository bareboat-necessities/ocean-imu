# BRMM stability proof handover

## Landing and starting point

The user authorized landing [PR #502](https://github.com/bareboat-necessities/ocean-imu/pull/502)
and [PR #504](https://github.com/bareboat-necessities/ocean-imu/pull/504) into main
on 2026-09-08. #502 merged as `14b4abc5b0ec0b24748c7acea9b3b8ac10558793`.
#504 was retargeted to that main revision; its merge-tree check was clean.
This handover is delivered by #504. Verify its merged status when resuming.
Both changes are separate from previously merged #500.

Start a new thread from freshly fetched main containing this document and
`tools/stability/ou3_brmm_p3_premises.py`. Read this handover and the research
ledger, then check CI on the landed revision. Do not start from an older
source/theorem or silently discard the stack's changes. Landing the research
work does not assert that every CI job or proof obligation has passed.

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

## CI at handover and outstanding validation

The first migration revision is `0db18980f4ff7e018e392a847ed4ba0d6ad24ef7`, tree
`de271057b08b870de9fb05f4099fc8518d7f3df3`. GitHub run 34172553221 passed the
source foundation, complete source, same-history execution and actual theorem
syntax build. Its Python quality job found one unused local, corrected in the
follow-up. A local reference hard-set paper-binding failure was also corrected
and the focused precondition test passed.

The corrected code revision is `1395e615107d26804f4fbe2ee08acdf1bd6e43ef`, tree
`b29fdb7c199456638141cb79a17340c3c36cf1c7`. At the landing checkpoint, run
34173154195 had passed source foundation and complete BRMM source; theorem and
same-history jobs were running, and P3 job 101897660299 was queued. **The new
BRMM P3 numerical rebuild is not reported complete at handover.** Follow the
landed revision's CI and record its actual H18/A21 delta, pivots and artifact
identity before declaring the conditional certificate rebuilt. Do not transfer
PASS from another source scope or a different code tree. The final landing
checkpoint changes documentation only.

Mandatory local `make all` fails because Eigen/Dense is unavailable. CI owns
the C++ build result. An independent OU validation gate reports stale replay
provenance for OU-II, OU-III, the tuner and the OU-III simulation driver. It
requires genuine full validation/robustness regeneration. The existing
`ou-full-evidence-branch.yml` workflow is opt-in and runs about 1150 simulator
replays; main has its automatic full-study path. Never repair this gate by
restamping hashes or bypassing the evidence contract. These are outstanding
validation obligations after the user-authorized landing, not proof of a
filter failure and not permission to promote P4/P5.
