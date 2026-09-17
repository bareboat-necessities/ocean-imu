# OU-III ALT current handover

rho has now been measured, and the declared joint24 contraction is falsified.
`ou3-alt-word-rho-feasibility.md` is the controlling result: an admitted ungauged
legal word -- one with no magnetic event -- carries an exactly invariant
unipotent `(theta_z, bg_z)` block, so its motion-block spectral radius is one and
no common storage with the declared `[e_ba, beta_true]` supply attains `rho<1`.
The obstruction is the formulation, not an enclosure, and no subordinate lemma
can move it.

The next task is therefore to choose a new formulation, not to finish the
remaining qualifications. `AGENTS.md` and `ou3-alt-proof-plan.md` are normative.
The independent P2/P3/P4/P5 route continues through
`ou3-brmm-main-handover.md`; its thresholds and gates are unchanged. The ALT
master still withholds storage search.

## Physical input theorem

Every admitted Live history now satisfies the actual binary32 API bounds
`max_i |gyro_body_i| <= 35 rad/s` and
`max_i |accel_body_i| <= 160 m/s²`, after calibration and body-axis mapping.
`ou3-alt-live-input-contract.md` gives the device evidence and exact profile.
Checks occur before execution and retain complete rounding cells. Startup
keeps the already selected BMI270/MPU6886 residual and temporal contracts;
these must not be selected again or silently widened.

The bounds enter the proof. The private Mahony integral has an all-prefix
rounding barrier 4096, corrected rate stays below 4132, and Euler norm sum
below 5000. Normalization returns the state to squared norm below 1.112,
giving vertical acceleration below 322 and the WPE envelope 512. Startup
memory is retained, with no finite capture deadline assumed. The normalization
bound includes the literal fused Newton correction on every finite norm cell.
State bounds do not imply bit identity between compiler modes or qualify
unused Euler-angle outputs, seed reachability, or the complete MEKF.

## Eleven qualifications

`tools/stability/ou3_alt_contraction/finite_master_guard.py` is authoritative.
Its inventory retains all eleven names; each value is derived from its
supplying proof. A closed library or supply row uses the named shared
execution premises. The independent firmware/compiler and source-reachability
gates prevent those premises from being silently discarded.

| Qualification | Current result | Remaining requirement |
| --- | --- | --- |
| 1. Timeout/alignment reachability | The 30,002-sample deadline is falsified on a source-audited native wrapper history; that history first reaches Live at 33,447 | A valid universal startup argument; later recovery of one history is insufficient |
| 2. Near-antiparallel Eigen SVD | Closed at the source-algorithm level: pinned 3x2 QR/Jacobi producer, both pivots/ranks/tiny tails, two-sweep termination, local same-tree FMA family, and first-seed bridge | Target whole-firmware/compiler/FCR correspondence remains a separate deployment row |
| 3. Universal startup reachability | Both commissioned profiles admit the missed-deadline history; explicit bad equilibria are available as conditional constructions | Reachability or avoidance of the bad startup set from actual construction |
| 4. Startup deployment supplies | MEMS all-prefix observer and configured frontend bounds, with actual state ancestry | Remaining startup initialization, Racc/MEKF arithmetic and target composition |
| 5. Target compiler/libm profile | Component compiler, Eigen, libm, scalar namespace, and WPE expression maps are pinned | Whole-firmware call graph and FCR.RM=0 initialization/preservation are not qualified |
| 6. WPE target selection | Closed for the pinned WPE machine graph: actual mixed FMA/MSUB sites map into the persistent FMA history | Whole-firmware compiler/FCR premise remains outside this component |
| 7. WPE uniform supplies | Closed for the declared initialized MEMS/dormant-guard envelope; frontend, candidate tau/sigma/R_S and positive powf range are bounded without a startup deadline | Startup-root and whole-firmware premises remain separate |
| 8. WPE log/exp libm | Closed for pinned newlib exp/log and IEEE sqrt error profiles under the named scalar/link premises | FCR/whole-firmware execution premise remains separate |
| 9. Qaxis exp libm | Closed for pinned exp error below 2^-24 and the Q-axis envelope | FCR/whole-firmware execution premise remains separate |
| 10. All-event arithmetic | Gradual underflow, input totality, and literal FMA state bounds are repaired | Source-uniform carried MEKF bias/state/covariance, solves and roundoff |
| 11. Complete 600-step word | Each mode carries its own CORE, watchdog, MAG and HOLD successors; actual rounded packets and goLive a_w reset are attached; the sealed ungauged startup CORE constructor is now represented | Universal startup root, gauged pending-yaw ancestry, quaternion lower shell, and rounded full CORE/Eigen totality |

## Startup evidence and its limits

`ou3-alt-startup-timeout-witness.md` supplies one analytic physical wave and
bounded yaw/residual history. Every native API packet is audited against that
same history, with zero BIAS0/1/2 driver, both commissioned sensor profiles,
no installed filter state, and no pre-Live magnetic calls. The guard remains
dormant. At both samples 30,002 and 30,602 the unchanged wrapper is not Live;
it subsequently recovers. This refutes the assumed deadline, not eventual
startup. Conditional inverted equilibria are not relabeled as reachable.

`ou3-alt-startup-direction-sampling.md` separately shows why the continuous
direction mean/primitive budget cannot be reused unchanged as a sampled budget.
Neither result permits a new input restriction for proof convenience.

## Retained finite-word structure

The source product retains corrected COMPLETE-BRMM, one BIAS0/1/2 history,
full 21-state covariance, joint24 motion/error/true-bias state, one Live/S
origin, frontend/tuner memory, asynchronous magnetic state and all H18/A21
edges. Four attitude charts cover every nonzero relative quaternion, including
south-heading entry, without a small-angle capture premise. The optional
wind-heel retarget and lever arm remain outside the declared scope.

Signed magnetic counters are safe by saturation on every finite prefix.
Certified empty startup magnetic histories, ungauged Live waiting and later
north acquisition are represented. External hold may preserve H18 indefinitely;
no eventual A21 transition is assumed. All source and machine predecessors
must belong to the same carried history.

The target arithmetic certificates intentionally retain the execution premise
`FCR.RM=0`. The pinned ISA leaves FCR reset undefined and the first FPU use
inherits the active register; absence of a startup writer is not an RNE proof.
The FCR audit is therefore recorded as an explicit deployment blocker rather
than silently promoted into the library/ compiler rows.

Read `ou3-alt-machine-core-continuation.md` for the remaining machine startup
root, and `ou3-alt-target-arithmetic.md` for exact target/library qualifications.
Read `ou3-alt-wpe-uniform-supplies.md` for the all-time WPE induction.
`ou3-alt-event-arithmetic-domain.md` retains the giant-gyro overflow as an
outside-domain regression; it is no longer the commissioned-input blocker.

The wrapper's 30,602-sample clock/error certificate remains a conditional
finite-horizon component. It cannot be used to prove that every admitted
startup lies in that horizon. Integer source ordinals, machine clocks and
literal branch predicates remain distinct. Finite clock lookup or an initialized
observer bound does not prove indefinite clock progress.

## Next controlling work

Pick a replacement formulation from the three routes in
`ou3-alt-word-rho-feasibility.md` -- gravity-quotient storage, a declared Normal
Live magnetic service class, or an independently bounded heading supply -- and
re-measure rho before any enclosure work. The gauged floors are already below
one with roughly `4.0e-3` distance per 3 s word, so for a gauging route the next
question is whether interval enclosure over 600 steps fits inside that margin;
measure the enclosure width before building it.

The eleven qualifications below stay open and correctly stated, but they are
subordinate to `rho_w`. Do not resume them as the main task: finishing all of
them would not reach the declared contraction. The machine Joseph update can
have a small positive covariance defect, and the literal first-order attitude
reset is not orthogonal; both charges must be retained by whatever formulation
replaces the falsified one. A complete word must still pass the finite-master
guard before any storage search.
