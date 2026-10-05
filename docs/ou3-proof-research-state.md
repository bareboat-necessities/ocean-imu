# OU-III proof research state

## Current hypothesis

The controlling objective is **all-time actually-applied MAGNETIC SERVICE for
the exact planar MOVING sensor-compatible pair**. Admission and exclusion are
both OPEN. Do not return to scalar-span bias separation or generic A21 entry
before resolving this question. Base inspected: main
`59959d1bf33a3c5f4f70331104546a6486468df6` (newer than PR #651).

The theorem remains unclosed: `theorem_closed=false`,
`regional_practical_stability_claimed=false`, `certified_capture_time=null`.
The manifest now has **24 PROVED, 4 CONDITIONAL, 1 OPEN**; the two new proved
items are local admission-cell algebra and planar guard inactivity in real
arithmetic. Neither is an all-time service certificate. Conditional obligations
remain `joint_kernel`, `physical_nominal_bridge`, `h18_release`, `compact_outer`.
The qualified local real-arithmetic result at sqrt(V)<=0.15 is retained, not
outer capture. Literal construction/capture/H18/refinement/release/A21 chronology
and physical S are inherited throughout.

## Evidence

**PROVED analytical / exact rational checks.** The planar pair has the same
complete sensor record for all time, with current MARINE/SLOW+FAST bounds and
constant tilt-induced BA difference. The same-prefix metric comparison is
max(V_plus,V_minus)>=1600(g sin(theta))^2>15.38649 on the specified SPD default
BA branch. It does not establish all-time magnetic admission.

The exact covariance parity remains 12+9. Service is reconstructed for the two
physical 2x2 +/- probes: after the orthogonal sum/difference transform each is
(I_E+I_O)/2, not min(lambda_min(I_E),lambda_min(I_O)). Actual acc and S updates
consume paired probe storage just as magnetic updates do.

`planar_service_cell.py` proves the same-P linked local prediction/S resolvent
identity, rank <=4/even and <=2/odd, and the fixed-target Frobenius
nonexpansiveness of the literal additive AW floor. AW synchronization is **not
Loewner monotone**. Actual reset is I+(1/2)[dtheta]x. Fixed-coefficient covariance
partial derivatives retain prediction, Joseph, reset and every active AW
marginal replacement; generated-mean/coefficient derivatives are separate open
ports. These identities enter the invariant admission cell, not a substitute
outer-tail theorem.

`planar_service_guard.py` proves all-time inactivity of the seeded default guard
for this exact source in real arithmetic: detector RMS <=0.000530491526<0.03.
It retains both detector stages, per-axis EW squares, slew and rail parking.
Float32 transfer remains open.

**FINITE DIAGNOSTIC ONLY.** A fresh 240-s full shipping replay gives minimum
actual transported/innovation-whitened one-second magnetic information
7.024764605642485, 3756 accepted magnetic updates and 102 disjoint service words.
The inherited final-20-s sweep checked 4001 sample roots with minimum
7.0275431864. Neither sweep covers all future histories/phases.

The new `planar_service_stream.py` exporter uses the true pre-prediction P and
literal F/Q/R_S, every correction H/R/innovation/K/PCt/residual, pre/post reset,
AW target/pending operation, and sample means/covariances/tuner outputs. An
untapped control matches every complete tail sample bit-for-bit. The 200--240-s
word has 8000 predictions/acc corrections, 1000 mag corrections, 294 S events,
and 381 pending AW synchronizations. Parity cross-blocks are exactly zero in
that replay. All 381 AW faces are active; observed target-minus-prior AW
minimum eigenvalue is 0.0008455758. This is not a uniform active-face certificate.

The corrected local S/prediction defects have maximum relative norms
0.0008993008/even and 0.0016318965/odd. Linked one-percent-cell formula evaluations
are 0.0016300211 and 0.0025795765. These floating local bounds are not complete
scheduler-word or outward-rounded bounds.

Whole 200--220-s covariance partial-derivative gains are 0.7426143369 and
0.8738212971 in the actual, different root/end covariance metrics, with all
AW derivatives retained. This is promising fixed-coefficient feasibility, not
a full joint derivative, periodic return or invariant-cell certificate.

**CONDITIONAL exact calculation.** `planar_service_mahony_tube.py` uses
z=(pitch error,10*integralFBy), G=[[10,-5],[-5,15]], effective feedback in
[.49,.51], and radius .006. Exact endpoint LDL checks give squared decrement
.00027; proposed forcing caps give charge 7841/12000000000, below radial margin
.00000081. Literal initialization and every-step float residual caps are not
bound yet. Raw quaternion norm is not assumed one.

## Joint phase-cell threshold calculation

A new fail-closed comparison calculation converts the existing literal 20-s
covariance partial derivative and exact conditional private-Mahony decrement
into a quantitative coupling target.  The worst parity covariance gain is
0.8738212970667966.  Ignoring additive charge only for the homogeneous
comparison, the exact .00027 one-step squared-norm decrement over 4000 samples
gives a 20-s Mahony homogeneous gain 0.582705763926933.  Thus a normalized
nonnegative two-block comparison [[a,b],[c,d]] is contractive whenever

    b*c < (1-a)(1-d) = 0.05265364544920159.

For symmetric coupling this is ||coupling|| < 0.22946382165649032.  This is a
FEASIBILITY THRESHOLD, not a bound on the generated mean/tuner coupling; those
ports are still open and the additive Mahony charge is retained separately in
the radius calculation.

The same calculation exposes the exact service perturbation budget.  Starting
from the carried physical +/- 2x2 floor 7.024764605642485 and mu_M=1, Weyl's
inequality leaves spectral-norm budget 6.024764605642485 for the complete
joint-cell information perturbation in each physical 2x2 block.  Therefore the
remaining every-placed-window proof does not need to reproduce the carried
floor tightly: it needs a same-history outward cell plus
||I_cell-I_carried||_2 < 6.024764605642485 uniformly.  Neither inequality is
promoted because generated coefficient couplings, continuous/root placement,
and forward containment are not yet certified.  See
planar_joint_phase_cell.py and planar-joint-phase-cell-threshold.json.

## Current limiter

The joint mean/covariance/frontend/tuner/clock cell is NOT constructed or proved
forward invariant. A 20-s source period is not an S or AW clock return. The cell
must include the planar mean, P_even/P_odd, private **raw** Mahony quaternion and
integral, guard, frequency/variance/period, staged/applied joint tuner tuple,
S elapsed, AW clock and pending target, and reference/refinement/gate state.
A tuple of independent observed marginal intervals is not a reachable family.

The inherited interval frontend was not source faithful: exact instead of fast
float inverse square root, omitted gravity subtraction, missing literal tilt
seeding/readiness, incorrect period/tuner ordering and unbound constants/clocks.
`planar_service_frontend_binding.py` gives an exact binary32 counterexample to
unit norm after normalization. The reference HistoryCell entry now fails closed
as **E_IMPLEMENTATION_FAILURE**, not a D subdivision failure. The old reference
recurrences remain available only as explicitly selected reference diagnostics.
The pre-normalization quaternion identity itself is unaffected.

## Failed approaches / DEAD_ENDS

* **E:** unpopulated prior-covariance tap made the old S-shift field zero, with
  hard-coded rather than actual R_S. Retired field is not evidence; use the new
  operation stream and linked shift audit.
* **E:** P<=Pupper does not imply ||PH'||<=||Pupper H'||. Exact SPD counterexample
  gives D11=500/3 against the old bound 1/20. Fixed with PSD block Cauchy and the
  linked rank factorization. This is not a physical counterexample.
* **E:** the orientation-free rank-three accelerometer outer product has false
  null directions. Use actual exported rows; its optional generic majorant is
  now an explicitly relaxed block-diagonal Young bound. S noise floors include
  the literal (.72,.50,1) axis factors.
* **E:** paired covariance/probe equality without measurement loss, an AW
  Loewner-monotonicity assumption, and an exact-normalization frontend replica
  are invalid. Corrected/scoped in code, status, manuscript and parity note.
* **D:** products of local relative-Frobenius AW norms lose whole-word
  cancellation. Per-AW odd log-product was 77.7623; one grouping refinement over
  complete between-sync words still gave 38.3687. Do not repeat this norm-product
  tactic; use the complete operator and preserve cross-blocks.
* **D:** generic-root scalar and anisotropic prediction-only product lower
  factors (3.759e-16 and 1.859e-13) cannot establish the service floor and omit
  intervening nonmagnetic losses if multiplied by a replay floor. Retired.
* **D:** scalar completed-square relaxation exceeded the absolute local-entry
  budget by more than six orders of magnitude. It is not shipping instability.
* The two-/three-epoch or longer-window physical BA separation mechanism is
  refuted by the exact constant-difference MOVING record. No interval sharpening
  repairs it. Do not remove physical BA forcing or strengthen excitation.

No genuine admitted shipping counterexample (A) has been established. The
absolute-entry impossibility remains conditional on full MOVING admission;
quiet absolute-entry obstruction is already retained independently.

## Retained facts

Preserve all manifest-proved operation/kernel/source lemmas, the full 21-state
implementation, literal covariance/gains/Joseph/reset/AW sync, BA/BG projection,
Mahony/reference acquisition, guard, frequency/variance, joint tuning and
progress-preserving S scheduling. Physical S is never reset by S=0 corrections.
The four conditional bridges have not been promoted. Quiet compatible-tube
boundedness, H18 release completion, quotient source-uniform action/supply,
every-prefix retention, nonlinear remainder, float32 and regime composition
remain open as recorded in the manifest/status.

## Alternatives

If all-time admission closes, the surviving physical gauge must remain in the
covariance-metric compatibility quotient/tube, with explicit gauge-to-transverse
C_Q alpha forcing. If literal service instead fails, prove that failure on the
same history before excluding the pair. Neither branch is selected yet.

## Next falsifiable experiment

Bind the conditional raw-Mahony pitch/integral cell and complete literal tuner
and clocks; do not use the retired ideal frontend output. Complete the joint
mean/coefficient sensitivity around the exported parity covariance operator,
including S-placement switches, AW active-face/target changes and reset effects.
Form an outward-rounded **joint** source/clock-phase cell and check containment
and every-prefix gates. Then sweep every placed one-second physical +/- service
root and prove its 2x2 floor exceeds one. A failure remains D unless it supplies
a fully admitted violating execution; no generic A21 detour is authorized.

## Verification scope

The source-pinned local OU-III suite passed all 691 tests. The regenerated
manifest/status/provenance validator passes with theorem_closed=false. The
240-s native export passed an untapped-control, bit-for-bit comparison of every
recorded tail sample and terminal state; the operation audit and covariance
partial-derivative verifier pass without theorem promotion. The authoritative
LaTeX study and appendix compile together with no overfull boxes.

The retired ideal frontend now causes the release and physical-subdivision
entry points to return an explicit E_IMPLEMENTATION_FAILURE, rather than
claiming a shipping HistoryCell. This is a deliberate fail-closed model-binding
check, not a proof of physical failure. The first local `make all` attempt failed
at missing Eigen headers. A retry uses the separately retrieved Eigen include
tree without changing repository build flags; the broad build is not counted as
a pass here. CI and broad-build completion are reported separately in PR status.

## Shipping-faithfulness handoff

1. Preserved: all shipping source files, constants, assumptions and quality
   gates; literal 21-state execution and inherited physical/filter history.
2. Relaxations: declared local covariance cells, fixed-coefficient covariance
   partial derivative, and a conditional frontend perturbation tube. None is
   promoted to a reachable shipping family.
3. Failures: E proof-model/instrumentation/inequality errors and D norm-product
   failures above; native network/source-access limitations are infrastructure E.
4. Admitted counterexample: none. MOVING admission/exclusion still OPEN.
5. Retained: qualified local theorem, exact pair/source/metric results and
   manifest-proved lemmas. Two new component lemmas do not close the theorem.
6. Next: actual joint causal phase-cell containment and physical 2x2 service
   floor, retaining mean, raw Mahony, both clocks, AW sync and every correction.
