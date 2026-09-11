# OU-III parallel ALT proof handover

## Resume from latest main in a new PR

Read `AGENTS.md`, this handover, `docs/ou3-alt-contraction.md`, the ALT section
of `docs/ou3-proof-research-state.md`, then `docs/ou3-brmm-main-handover.md`.
The original P2/P3/P4/P5 route remains independently continuable. Do not change
its implementation, proof tools, gates or shared source contracts to help ALT.
ALT does not consume its PASS labels as a certificate of the ALT master.

## Target and immutable contracts

Keep the joint true-minus-estimate state
`z=(e_theta,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta_true)` and coercive storage with
motion/bias cross terms. Seek
`V_next <= rho V + w^T Gamma w + c^T Beta c`, `0<rho<1`, with a quantitatively
useful bound. Only independently bounded physical/reference/bias quantities
may enter c. Unknown motion errors are not bounded input ports. Strict
homogeneous contraction of all 24 coordinates is ruled out by held/constant
bias directions; this does not rule out the desired practical ISS theorem.

Preserve corrected COMPLETE-BRMM: acceleration cap 8.8 m/s^2, centered primitive
cap 1100 m*s, one-time Live S origin, no position reanchor, zero lever arm and
transparent dormant vibration guard. Keep actual anisotropic R_S, same-signal
WPE -> sigma -> tau -> T_S -> R_S, staged tuner/scheduler state, H18 and A21,
actual Joseph/reset/projection, separate BIAS0/BIAS1/BIAS2 ancestry, and the
finite-precision obligation. P3 stays at 1e-18. No replay fitting, finite-seed
source qualification, artificial domain reduction, or covariance consistency
as hard entry membership.

## What is executable now

`core.py` retains inverse-free measurement/increment algebra, the shared-driver
joint24 bias lift, projection IQC, sign-correct master assembly and structural
anti-shortcut tests. The old frozen-map observer remains a separate baseline.
It does not bind endogenous gain dependence or the actual full joint24 word.

The new actual-word experiment is:

```
python tools/stability/ou3_alt_contraction/shipping_probe.py --output /tmp/alt --family 0
python tools/stability/ou3_alt_contraction/shipping_probe.py --output /tmp/alt --family 1 --skip-build
python tools/stability/ou3_alt_contraction/shipping_probe.py --output /tmp/alt --family 2 --skip-build --detailed-direction -5
python tools/stability/ou3_alt_contraction/word_diagnostic.py --directory /tmp/alt/BIAS2 --output /tmp/alt/BIAS2/diagnostic.json --solves
```

The probe uses a temporary, observation-only header overlay. Removing its
callback lines must reproduce the original header exactly. It forks the entire
running process at literal word entry, preserving private covariance, frontend,
tuner and guards. Initial error directions include held bias; physical-bias
root variants restart from boot, never from a new root at a word boundary.
The analytic source is deterministic, not a captured wave replay. Its all-time
kinematic/bias bounds do not qualify the complete BRMM/PE/startup domain.

There are 225 actual three-second tapes: three bias families, three natural
H18/A21/release checkpoints, and 25 executions per checkpoint (baseline,
21 error perturbations, three boot-root variants). The analyzer retains all
24 error/truth coordinates. It does not assemble these finite secants into a
Jacobian, or treat unrelated boot-root covariance histories as its last three
columns. Divergent event words are reported, not silently aligned.

Detailed tapes bind actual N/S/r/K/P and literal pre/post means to the finite
solve identities. The analyzer separates gain-solve defects from float
multiply/add defects. High-precision evaluation at 80/120 digits checks these
identities and the endpoint chart on recorded binary32 operands; it does NOT
re-execute the nonlinear shipping algorithm at arbitrary precision or enclose
deployment roundoff. Run the 23 ALT tests with
`python -m unittest discover -s tests/ou3_alt_contraction -v`.

## Evidence and failed candidate

The predeclared diagnostic metric uses block scales
`(.1,.01,1,1,10,1,.1,.1)` and normalized bias/truth cross term 1/4. It was not
fitted to the tapes. Its worst tested *unsupplied finite-secant* ratios are:

| Family | H18 | A21 | H18 -> A21 |
| --- | ---: | ---: | ---: |
| BIAS0 | 21.5559 | 9.96772 | 9.74203 |
| BIAS1 | 21.3234 | 9.92484 | 9.73857 |
| BIAS2 | 21.5756 | 9.93982 | 9.73054 |

The limiting directions are e_S,y, e_theta,x and e_v,z respectively. This
rejects that metric as an unsupplied incremental contraction witness, not the
practical ISS theorem or every common joint storage. No justified supply has
been subtracted. Ratios are not comparable to the old observer's covariance-
weighted frozen-map ratios because the metric and state experiment differ.

Nonzero nominal residual, delta_N*q0, delta_S*q0 and covariance increments are
now measured directly. The H18 attitude-direction probe has an omitted
frozen-gain correction norm about 6.8281e-6 in raw mixed coordinates. This is
not a uniform performance bound. The recorded-gain/ideal-solve discrepancy
reaches about 8.2452e-9, while the full literal mean-update defect reaches
9.3577e-7 after multiply/add rounding is included. These raw mixed-coordinate
maxima are measured defects, not certified supplies. The uninstrumented control
reproduces all 50,775 recorded wrapper states across 75 BIAS2 words bitwise
under the same compiler/options; this is not deployment qualification.

See the current ALT research ledger for the failure classification, critic pass
and alternatives. Do not tune the rejected metric to these tapes or begin
interval subdivision from this evidence.

## Next falsifiable experiment

Complete the missing **source-uniform graph attachment**, using the executable
probe as a correspondence check, not a universal source certificate. Specify
which covariance/frontend/guard increments are initial internal states and
which are generated by the same history. Prove their consecutive-word
compatibility; do not silently set them to zero again at each word.

Before implementing the next master, explicitly choose between (1) a common
joint-storage master with analytically constrained internal increment ports,
(2) storage augmented by the necessary internal/covariance increments, and
(3) a same-trajectory physical-error identity retaining reference forcing and
legitimately avoiding comparison of independent Riccati histories. The chosen
formulation must account for the measured limiting quantities. Then execute a
high-precision complete-word feasibility diagnostic of that formulation over
an analytic source family, without fitting to the retained tapes. Only a
feasible master justifies rigorous outward/interval work.

## Still open, fail closed

The shared aggregate bias-supply builder still fails fresh-Live hard membership
with its original coordinate/cover errors. Individual BIAS0/1/2 definition
validators pass, but startup and hardware admission do not follow. Do not
bypass the aggregate failure. The original PE/vector-domain issue is separate.

Uniform compatible storage, justified supply coefficients, source cover,
nonlinear/hybrid graph closure, every-prefix chart retention, useful ultimate
bound and basin, actual Mahony/proxy capture on both measured-period and
prior-frequency-timeout paths, release landing, and deployment finite precision
remain open. Observed natural entry times are not universal capture bounds.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Existing `P4_PASS=false` and `P5_MAY_START=false` are unchanged. Green experiment
CI means that diagnostics executed, not that either theorem route closed.
