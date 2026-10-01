# A21 staged entry and prefix retention after field-axis exclusion

Source: PR #637 after the radius-local field-axis exclusion. No estimator,
assumption, tuning or quality-gate change is proposed.

## 1. Direct H18 -> V<=.15^2 entry is not a valid theorem target

While accelerometer-bias updates are held, shipping explicitly zeros every BA
cross covariance with the base and LIN blocks. The BA covariance is seeded
with sigma_bacc0=.004 m/s^2. Enabling A21 only applies

    P_ba,ii <- max(P_ba,ii, sigma_bacc0^2),

so there is no release covariance inflation capable of hiding an arbitrary
held physical BA error.

The admitted physical residual satisfies

    ||b_a|| <= B_a = .22516660498395405 m/s^2.

For an admissible held execution with b_hat_a=0, the release storage contains
the exact decoupled BA contribution

    V_ba = ||b_a||^2/.004^2.

At the amplitude envelope this is >3168, whereas the field-axis local ball is

    r_FA^2=.0225.

Thus the present assumptions cannot imply universal direct release into
V<=.15^2. This is not a filter failure and does not refute eventual A21
capture. It proves that "release already lies in the final local ball" is the
wrong entry lemma.

## 2. Homogeneous every-prefix retention inside the local ball is closed

For the covariance-matched homogeneous comparison, every literal prediction

    P^- = F P^+ F' + Q, Q>=0

is storage-nonexpansive, every applied Joseph correction satisfies the exact
measurement-energy identity, and a reset applied congruently to error and
covariance preserves storage. PSD covariance-floor/release events can only
decrease the inverse-covariance storage of a fixed error. Therefore at every
regular real-arithmetic prefix

    V_prefix <= V_root.                                     (ER1)

Consequently V_root<=.15^2 implies homogeneous every-prefix retention in the
same ball. The BA estimate projection is inactive there because the much
larger already-proved projection guard is sqrt(V)<=6.

ER1 does not include finite nonlinear injection mismatch, physical/source
mismatch or arithmetic residuals. Those enter the prefix supply below.

## 3. Correct staged entry coordinate

At A21 release BA is decoupled by the preceding held operation. Partition

    e=(e_o,e_ba),   P_release=diag(P_o,P_ba)

at that instant and define

    V_o=e_o' P_o^-1 e_o.

The field-axis invariant-set proof depends on e_aw in e_o and not on BA:
the accelerometer attitude Jacobian uses the nominal CoG vector a_hat_w-g,
while BA is a separate measurement column. Therefore the LaSalle
field-alignment exclusion may be applied to a candidate OUTER/BA-quotient
storage provided the source-audited AW marginal implication

    ||e_aw|| < 4.06 sqrt(V_o)                               (ER2)

is retained for the corresponding Schur/decoupled covariance. At the literal
release boundary it is immediate from block diagonality. After release, BA cross covariance regrows. Eliminating BA from the STORAGE
quadratic means minimizing over e_ba:

    V_o,elim(e_o) := min_z [e_o;z]' P^-1 [e_o;z]
                   = e_o' P_oo^-1 e_o.                     (ER3)

Thus the eliminated covariance is the OUTER MARGINAL P_oo. The covariance
Schur complement P_oo-P_ob P_bb^-1 P_bo instead corresponds to conditioning
on/fixing the BA coordinate and is NOT equivalent to eliminating BA from the
quadratic. At the literal held-to-active release boundary P_ob=0, so the two
happen to coincide there. After release they must not be conflated.

This eliminated storage is the correct quantity: a large released BA error
must not destroy already-local attitude/AW geometry merely because full V is
large.

The next staged theorem is therefore:

(A) prove V_o,elim<=.15^2 entry/retention through capture/H18/release;
(B) on that quotient-local A21 class, use the field-axis LaSalle exclusion to
    obtain finite-window strict quotient dissipation;
(C) carry the bounded BA error through its literal active OU prediction,
    accelerometer correction and dissipative estimate projection until full
    V enters .15^2;
(D) thereafter use the full local theorem.

This is not a new architecture; it is block elimination of the held coordinate
already required by the H18 complement proof.

## 4. Nonlinear/source prefix and endpoint closure

Let q=sqrt(1-eta_D) <1 be the homogeneous finite-window norm factor supplied
by compactness on the quotient/full local class. For the literal finite-error
composition define, on radius r,

    E_W(r,d) = complete-window additive sqrt-storage supply,
    G_k(r,d) = additive supply accumulated to prefix k,
    g_k(r)   = homogeneous prefix norm gain.

The exact operation composition has

    sqrt(V_next) <= q sqrt(V_root)+E_W(r,d),                (ER4)
    sqrt(V_k)    <= g_k(r) sqrt(V_root)+G_k(r,d).           (ER5)

For the homogeneous covariance-matched comparison g_k<=1 by ER1. Finite
nonlinear coordinate/reset remainder may instead be kept in G_k; it must not
be hidden in g_k without a proved Lipschitz coefficient.

A full local self-map at r_FA=.15 therefore requires

    E_W(.15,d) <= (1-q).15,                                 (ER6)
    G_k(.15,d) <= (1-g_k).15 for every prefix k.            (ER7)

In the ER1 normalization g_k=1, ER7 shows why a zero-slack use of the same
.15 boundary is inappropriate for nonzero additive supply. Choose an inner
root radius r_in<.15 and require

    g_k r_in+G_k(.15,d) <= .15 for every k,                 (ER8)
    q r_in+E_W(.15,d) <= r_in.                              (ER9)

Equivalently

    sup_k G_k <= .15-r_in,
    E_W <= (1-q) r_in.                                     (ER10)

These are the exact retained-radius inequalities to be populated by the
existing actual-gain finite-error operation supplies. They compose sensor/model
defects, curvature, reset remainder, projection and arithmetic once, on the
same history.

## 5. What is and is not closed

CLOSED:
- the field-axis zero-dissipation candidate is excluded for the outer/AW local
  geometry at r=.15;
- homogeneous every-prefix storage retention once inside the full local ball;
- direct H18-to-full-ball entry is analytically ruled out as a valid universal
  target by the literal held BA covariance/error scale;
- the correct inner/outer prefix fixed-point inequalities ER8--ER10.

OPEN:
- quotient-local entry V_o,elim<=.15^2 at release;
- source-uniform quotient compactness/dissipation after BA elimination;
- a quantitative eta_D (hence q);
- exact E_W and G_k bounds small enough for ER10;
- finite active-A21 time until the BA coordinate makes full V<=r_in^2;
- float32 totality.

The next calculation should attack (A), not retry direct full-V release:
evaluate the actual capture/H18/release OUTER Schur storage with BA removed.
If that outer release set fits .15, the field-axis result supplies the missing
quotient invariant-set exclusion while the active BA transient is handled as
a bounded internal coordinate. If it does not fit, identify the outer
component consuming the radius before changing assumptions.


## 6. Literal carried release audit: attitude already excludes .15 entry

The dedicated read-only diagnostic
`release_outer_storage_diagnostic.py` snapshots the unchanged shipping
construction history at the first sample for which
`acc_bias_updates_enabled()` is true. It is finite carried evidence, not a
source-uniform release theorem.

Observed literal release:

    active_step = 36008
    t_release = 180.0399959757924 s
    tilt = 0.14215579822992427 rad = 8.14492727188796 deg

The release attitude covariance marginal is

    [[8.08563254395267e-6, 0, 4.691878984885989e-6],
     [0, 3.535204314175644e-6, 0],
     [4.691878984885989e-6, 0, 8.413949217356276e-6]]

with eigenvalues

    3.5352043141756427e-6,
    3.555041007836594e-6,
    1.2944540753472348e-5.

At release held BA is decoupled, so BA elimination and the outer release block
coincide. To identify a non-BA obstruction without needing physical v/p/S
coordinate reconstruction, minimize the outer quadratic further over every
non-attitude outer coordinate. The exact covariance identity leaves

    V_outer,elim >= theta' P_theta,theta^-1 theta
                  = 5716.295063726523,

hence

    sqrt(V_outer,elim) >= 75.60618403098071.                (ER11)

This is overwhelmingly outside r_FA=.15 (r_FA^2=.0225). The radius consumer is
therefore already ATTITUDE; BA elimination does not repair direct release
entry on this carried history.

A useful scale consequence follows without using the observed attitude
direction. Since lambda_max(P_theta)=1.2944540753472348e-5, membership in
V<=.15^2 would NECESSARILY require

    ||theta|| <= .15 sqrt(lambda_max(P_theta))
              = 5.3968e-4 rad ~= .03092 deg                 (ER12)

at this covariance. The observed 8.145 deg is over two orders of magnitude
larger in angle.

ER11 is not a source-uniform counterexample to eventual capture: this stress
history is already documented as not certifying all-time MAGNETIC SERVICE and
does not satisfy the desired captured release premise. It DOES answer the
literal carried-release diagnostic and shows that the next theorem cannot be
"stage flags imply V_outer<=.15^2". The missing entry mechanism must establish
subsequent A21 attitude/storage contraction from a substantially larger
release set until the trajectory enters the .15 local LaSalle ball.

Thus the proof needs two nested A21 regions:

1. an OUTER A21 capture/retention region large enough to contain certified
   release attitude/storage;
2. the INNER r_FA=.15 region, where field-axis invariant-set exclusion gives
   strict local dissipation.

The next analytical obligation is to use the zero-dissipation exclusion on
the inner boundary together with dissipativity on the compact annulus
.15 <= sqrt(V_outer,elim) <= R_outer to prove finite entrance, or to derive a
larger field-alignment exclusion radius valid on that annulus. Directly
asserting .15 entry at release is withdrawn.


## 8. Universal outer -> inner full-state entry is false under the present contract

The shaped-storage bootstrap cannot prove universal finite entrance into the point inner ball on all certified A21 continuations. The obstruction is the source-audited stationary attitude/accelerometer-bias gauge, now applied directly to the proposed entry theorem.

At physical rest the noiseless sensor map is

    gyro = b_g,
    f = -g Q e_z + b_a,
    m = Q B.

A constant-rest member suffices. Take Q_alpha to be a constant rotation by alpha=1/1000 rad about B, b_g=0 and b_a=g(Q_alpha e_z-e_z). Then p=v=a=0, all bias rates vanish, ||b_a||=2g sin(alpha/2)<.009807 m/s2 << B_a, and the accelerometer, gyro and magnetometer packets are exactly the quiet nominal packets forever. Complete stillness may persist arbitrarily long, so no T_E/theta_E excitation is owed. The actually applied quiet magnetic-service floor remains above the required value. The existing quiet covariance comparison gives P_ba,ba <= I/1600; hence on this distinct physical truth

    sqrt(V) >= ||e_ba||/.025 > .392 > .15.

Because the shipping estimator, tuner, covariance, magnetic reference, scheduler and gates are deterministic causal functions of the common packet history, the nominal and rotated physical truths generate the same complete filter history. Consequently no storage constructed solely from that history can force both physical errors to the point set V<=.15^2. OU forgetting, S=0 pseudo-updates, AW covariance synchronization and active BA prediction do not add physical information: the BA error recurrence contains the physical supply (1-phi_b)b_a+w_a, and magnetic corrections are blind to rotation about the magnetic axis. A signed/windowed AW primitive cannot remove this stationary gauge.

Therefore the requested implication

    every certified H18/A21 release -> eventual retained V<=.15^2 entry

is FALSE under the current assumptions. This is an identifiability/theorem-target failure, not a failed scalar bound and not evidence for changing runtime behavior. Do not interval-refine a shaped point storage against this family.

The controlling nested proof remains one architecture but must be regime/quotient aware:

1. define W_out on distance to the measurement-compatible attitude/BA class, with BA eliminated by the outer marginal rather than conditioned;
2. prove stationary A21 retention/practical stability to that class while carrying the literal tuner/covariance/service history;
3. on a complete admitted moving episode with the declared attitude-span excitation, prove that moving information collapses the stationary gauge and gives finite retained entry into V<=r_in^2, r_in<.15;
4. then invoke the existing radius-local LaSalle theorem and finite-error retention inequalities.

This does not strengthen MARINE MOTION. If an execution remains at complete physical rest forever, its theorem conclusion is practical stability to the stationary consistency class, not convergence of an unobservable physical attitude/BA representative to zero covariance-metric error.

Next falsifiable calculation: construct the exact stationary gauge tangent

    K_stat = {(delta theta,delta b_a): delta theta || Q B,
              delta b_a = g delta(Q e_z)}

at the carried A21 root, project the complete literal corrected-word action to K_stat^perp, and test whether one complete same-history T_E moving window supplies strictly positive action on the remaining gauge coordinate using the literal coupled tau/sigma_aw/R_S/T_S chronology. If the high-precision margin is not comfortably positive, identify the surviving same-history gauge before any interval refinement.
