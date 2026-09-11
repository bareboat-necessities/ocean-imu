# OU-III parallel ALT proof handover

## Where to resume

Continue ALT in a NEW PR from latest main. Read, in order: `AGENTS.md`, this
file, `docs/ou3-alt-contraction.md`, the ALT section of
`docs/ou3-proof-research-state.md`, then `docs/ou3-brmm-main-handover.md`.

The existing P2/P3/P4/P5 proof is an independently continuable route. Do not
rewrite it, make it depend on ALT, or consume its PASS labels as an ALT master
certificate. Shared source/runtime facts require their actual hypotheses and
same-history ancestry. The shipping implementation and shared proof machinery
are unchanged by the ALT attachment experiment.

## Immutable contracts and target

Keep corrected COMPLETE-BRMM, including `||a_wave|| <= 8.8 m/s^2`,
`D_S <= 1100 m*s`, one-time Live S origin, zero lever arm, actual anisotropic
R_S, same-signal WPE -> sigma -> tau -> T_S -> R_S, staged tuner/scheduler,
Joseph/reset/projection, BIAS0/BIAS1/BIAS2, H18 and A21. P3 stays `1e-18`.
No replay fitting, finite-seed qualification, wordwise S reset, covariance
consistency as hard entry, operating-domain shrinkage, or filter retuning.

Use true-minus-estimate joint24

`z=(e_theta,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta_true)`.

The target is coercive joint storage with motion/bias cross terms and justified
supply, not strict homogeneous contraction of persistent bias coordinates:

`m||z||^2 <= V(z,xi) <= Mbar||z||^2`,

`V_next <= rho V + w^T Gamma w + c^T Beta c`, `0<rho<1`.

Only independently bounded physical/bias/reference quantities belong in c.
Unknown motion or covariance-memory errors cannot be relabeled as bounded
inputs. A finite but enormous ultimate bound is not a useful Live basin.

## Executable attachment

`core.py` retains the inverse-free/IQC algebra, joint bias prediction and
projection sector. `shipping_graph.py` now derives conditional N/S/residual
relations from the same pre-state for accepted accelerometer, magnetometer,
and S updates. It does not infer H from P^-1 N. H18 masks gain bias columns/rows
while retaining bias uncertainty and cross terms in S. Actual R_S is read,
not replaced by a scalar nominal value.

`shipping_word.py` and `tests/ou3_alt_contraction/shipping_word.cpp` execute the
full default Fusion runtime with a continuous analytic physical source and
separate bias examples. Actual startup supplies baseline H18, H18->A21 and A21
roots. Host fork preserves all private runtime memory at each root. Tilt,
held-bias, and covariance-only finite perturbations yield 27 paired 600-step
words, with nine baseline words. The injected roots are NOT proved reachable.
The analytic examples satisfy checked kinematic/bias parameter bounds, NOT
full COMPLETE-BRMM/PE or hard-entry admission.

Temporary headers add observation hooks only; repository runtime files are
not edited. The optional uninstrumented control requires bit-identical
sample-boundary records for the same compiler/executions. Capture includes
actual P/N/S/K/residual, gain masks, mean update, Joseph, quaternion/covariance
reset, projection, attempts/returns, prediction/floor and bias-mode events.
Private frontend evolution is executed, not frozen. Emitted frontend summaries
are not a complete analytic recurrence or source-uniform reachable cover.
The default continuous magnetometer hard-iron correction and learned magnetic
reference are retained in the physical-input binding.

The shipping implementation solves each gain row before multiplying by r.
For its LDLT Lower view, write

`S_i^L K_i^T = N_i^T + E_i`,

`S_1^L delta_K^T + delta_S^L K_0^T = delta_N^T + delta_E`,

`delta(Kr) = K_0 delta_r + delta_K r_1`.

E is a measured solve defect, not a uniform roundoff enclosure. Raw S is kept
separately for the literal upper-triangle Joseph update. The driver
`w_b=beta_next-phi_true beta` is used once in both truth and bias-error graphs;
`(phi_true-phi_hat)beta` is retained. Conditional identities are useful
attachment primitives, not a complete source-uniform word.

Run on a Linux host with Eigen, NumPy, SciPy, mpmath and a C++20 compiler:

```sh
python -m unittest discover -s tests/ou3_alt_contraction -v
CXX=clang++ python tools/stability/ou3_alt_contraction/shipping_word.py \
  --work /tmp/ou3-alt-finite-word \
  --output /tmp/ou3-alt-finite-word.json --verify-uninstrumented
```

Use `--eigen PATH` when Eigen is elsewhere. The independent ALT workflow
retains the report, compiler/source fingerprints and manifest; full traces
are reproducible in the work directory. This experiment needs no downloaded
wave dataset. The older frozen-map diagnostic remains a separate baseline.

## Result and precise blocker

All three bias definition validators and the conditional operand regressions
are exercised. The covariance-only probes start with delta-z24=0 and end with
nonzero state differences on the same physical history. Thus an incremental
master cannot allow arbitrary covariance-memory differences while measuring
only delta-z24 and giving no memory supply. A nonzero nominal residual makes
the omitted endogenous gain term material even when delta-r initially vanishes.

This is a dependency/representation obstruction to that relaxed graph, NOT a
physical counterexample: the perturbed covariance roots are not admitted by a
source-uniform reachability theorem. It was anticipated by the warning that xi
must not be frozen. Do not repeat this probe with smaller intervals as a new
proof attempt. A source-dependent metric alone cannot make V(delta-z=0)
positive for an omitted memory difference.

No source-uniform rho, maximizing generalized direction, per-operation storage
margin, useful ultimate bound, or basin has been obtained. They remain null,
not manufactured from mixed-unit Euclidean norms. The 80/120-digit checks are
of recorded finite row-solve identities, NOT high-precision nonlinear word
re-execution. The old frozen-map ratios likewise do not supply that master.

## Next falsifiable experiment

1. Keep the new actual operand bindings and choose a closed master before any
   interval refinement. Compare three distinct routes: memory-augmented
   incremental storage; same-trajectory error-to-truth dissipativity; or a
   source-reachable pair invariant tying covariance/frontend increments to
   state and physical inputs. The current priority is a same-trajectory
   error-to-truth graph, which does not compare arbitrary covariance histories.
   Keep common joint24 storage plus bounded neutral supply as the first metric
   search; do not silently jump to independently fitted wordwise metrics.
2. Derive that graph's physical reference defect, especially S=0 forcing,
   actual covariance evolution, nonlinear attitude and frontend/guard ancestry.
   Truth is not an estimator solution. Neither physical forcing nor an omitted
   memory coordinate may be charged to an unproved bounded error port.
3. Once the complete master is closed, run its high-precision H18/A21/edge
   feasibility test with maximizing directions, supply coefficients, signed
   operation margins and distance to rho=1. Current finite probes do not meet
   this requirement. Obey the two-strike rule before refining a failed method.
4. Only if feasible: outward source-uniform verification, common/compatible
   coercivity, every-prefix chart retention, finite precision, useful basin and
   ultimate bound. Separately close actual Mahony/proxy capture for measured
   period takeover and prior timeout, plus H18->A21 landing.

Fresh-Live aggregate hard membership remains open; individual bias definitions
do not close it. The existing PE/vector consistency blocker stays on the other
route unless a shared lemma is actually required.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Existing `P4_PASS=false` / `P5_MAY_START=false` remain untouched. Green ALT CI
means the experiment and its regressions ran, not that stability is proved.
