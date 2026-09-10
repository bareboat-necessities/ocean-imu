# OU-III shipping Live-entry graph

## Scope and status

The target remains finite-time capture by the actual startup implementation,
followed by regional practical ISS of the 18 motion-error states and bounded
accelerometer-bias error for every admitted COMPLETE-BRMM/BIAS0/BIAS1/BIAS2
continuation. The entry identities below do not establish that theorem.
P4/P5 are not promoted. P3 delta remains 1e-18. No shipping arithmetic,
filter, source definition, R_S semantics, or quality gate is changed.

The distinction between a **fresh outer wrapper** and an arbitrary call to
an already-driven inner filter is essential. The executable paths are in
`src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h` and
`src/kalman_ou_iii/Kalman3D_Wave_OU_III.h`.

## Executable initialization and transition audit

The core constructor zeros the entire nominal extended state and initializes
P separately. The normal outer wrapper calls the frontend-only update during
Bootstrap, and calls the core-driving update after Live. Thus the fresh
wrapper's first Live handoff inherits zero v, p, S, a_w, gyro-bias and
accelerometer-bias estimates from construction. It does not integrate the
core navigation states while learning the proxy and tuner.

`goLive()` initializes attitude from the supplied proxy quaternion and enters
Live. It does **not** zero v, p, S or a_w. Consequently a generic `goLive()`
call after prior core propagation preserves those previously propagated
values. The constructor-zero argument cannot be used for that path, for a
subsequent Live word boundary, for a tilt reset, or for A21 release.

| Transition/call | Nominal navigation and bias means | Covariance and auxiliary state |
| --- | --- | --- |
| Fresh core construction | All extended means zero; quaternion identity | Block-seeded P; initial v/p/S standard deviations 1 m/s, 20 m, 50 m*s; scheduler elapsed zero |
| Bootstrap frontend-only updates | Core means do not propagate | Mahony/proxy, period estimator, bandpass, sigma statistic, tuner and related frontend state evolve |
| Cold/TunerWarm/TunerReady | No blanket nominal-state reset at these stage changes | Stage/counter and tuning operations must be retained; TunerWarm is not the H18 Live mode |
| `initialize_from_attitude()` | Attitude replaced; nominal attitude-correction coordinates zero; v/p/S/a_w/b_g/b_a preserved | Attitude covariance reseeded; attitude/bias crosses cleared; base/linear crosses cleared by the shared helper; no full P restart |
| First fresh-wrapper `goLive()` | Zero linear means inherited, not re-zeroed; proxy attitude installed | Live tuning applied; a_w covariance reseeded and its cross-covariances cleared; learned frontend/tuner state retained |
| H18 to A21 bias release | Navigation, attitude and bias means preserved by bias enable | Enabling b_a updates floors its diagonal covariance; does not reinitialize the whole P or the linear states |
| Live tilt reset | Attitude reinitialized preserving yaw; navigation and bias means preserved | Attitude-related covariance reset operations apply; this is not a constructor-equivalent restart |

`set_initial_linear_uncertainty()` is a covariance operation, not a nominal
state reset. `reset_aw_covariance_to_stationary()` clears a_w covariance
crosses and replaces that marginal; it does not zero the a_w estimate.
Disabling accelerometer-bias updates clears its cross-covariances while
preserving the bias mean. Enabling updates applies the configured marginal
variance floor. The actual scheduler retargeting operation must be carried
with tuning changes; it cannot be replaced by a newly chosen phase.

The outer magnetometer initialization/refinement and external bias-hold
release are additional hybrid operations. They must not be silently
identified with the inner stage enum. Likewise the timeout handoff path
cannot be treated as a guard explicitly requiring a usable period: the
quality and timeout branches have different conditions. Any assertion that
all reachable timeout handoffs already satisfy the period/heading premises
requires its own reachable-state proof.

## Answers about v, p and S

For a fresh normal wrapper, at the instant of its first Live handoff:

1. Estimated v is zero, inherited from construction.
2. Estimated p is zero, inherited from construction.
3. Estimated S is zero, inherited from construction.
4. The common nominal reset is core construction, not Live entry or A21.
5. Frontend learning does not propagate the core states on this path. A
   previously driven inner filter is a different entry family.
6. P is not generally restarted at Live: the attitude, a_w and bias-mode
   operations described above change selected blocks/crosses. Covariance
   magnitudes are not hard true-error bounds.
7. The shipping handoff does not receive physical truth and does not perform
   a physical truth-coordinate rebase. The physical/reference gauge must
   therefore be established by the admitted source definition, not inferred
   from assignments to estimated state.
8. In the unchanged physical gauge, e_v=v_true(t_L), e_p=p_true(t_L),
   e_S=S_true(t_L), and e_aw=a_true(t_L) on this fresh entry path.

Neither p_true(t_L)=0 nor S_true(t_L)=0 follows just from construction.
No covariance-confidence ellipsoid supplies the missing physical admission.

## The exact correlated entry constraint

Let x=(v,p,S,a_w), e=x_true-x_hat, and retain the *same* physical root on
both sides of the augmented state. At fresh entry:

    [e; x_true] = G x_true,       G = [I; I].

If H selects S, the literal S=0 innovation is

    r_S = -S_hat = H e - H x_true = [H,-H] G x_true = 0.

This is an exact identity. It is not a small-error approximation and does
not require S_true itself to vanish. An independent e_S ball detached from
S_true admits innovations that this entry path cannot produce. Conversely,
a physical source with nonzero S_true can produce a nonzero e_S at entry;
removing all S-error energy would also be wrong.

For a same-history finite gain K, the additive **linear-state rows** of an
S=0 update have augmented map

    D(K) = [[I-KH, KH], [0, I]],
    D(K)G = G.

The identity holds for any substituted gain, hence for the gain derived
from the actual same-history P/H/R. It does not permit independent K boxes
in a contraction certificate. Quaternion injection/reset, bias projection
and other sensor innovations are not replaced by this linear-row identity.

For a shipping linear prediction F, define

    d = x_true_next - F*x_true

from the same admitted physical continuation. Then

    [e_next; x_true_next] = diag(F,F)[e; x_true] + G*d,
    diag(F,F)G = G*F.

The same d must occur in both rows. Treating these copies as independent
forcings destroys the cancellation. Binary32 defects require explicit
additive arithmetic channels; the exact-arithmetic identity is not itself
a finite-precision enclosure.

For an augmented quadratic storage M, the entry restriction is the exact
pullback G^T M G. A residual-energy quadratic vanishes on the entry graph,
but true-error energy generally does not. The physical root still requires
its COMPLETE-BRMM admission constraints. The implementation in
`tools/stability/ou3_p4_correlated_entry_graph.py` checks exact polynomial
coefficients of these identities using rational arithmetic. Its finite
basis checks exhaust coefficients of affine matrix polynomials; they are
not sampled physical-history certification.

## The 300 m*s factor and the frozen diagnostic

The independent 300 m*s entrance factor is not established by a constructor
standard deviation of 50 m*s. Such a numerical relationship, even if present,
is not a provenance argument or a deterministic true-error guarantee.
The original introduction of 300 and its physical admission still require
explicit history/source evidence. No replacement radius is chosen here.

The frozen H18 envelope near 4.5788 versus chart bound 1 shows failure of
that frozen-map certificate on its chosen independent entrance factors.
It does not prove that an admissible nonlinear shipping trajectory exits
the chart: extrapolation outside the chart is not a valid nonlinear
counterexample. Likewise the lower frozen numbers after omitting the
independent S factor are not a source-uniform P4 certificate.

The corrected question is whether the full entry/source graph, followed
through every actual prediction, sensor update, adaptation/commit and
Joseph/reset operation, satisfies the retained-domain inequalities. Zero
innovation at entry must not be imposed again after other measurements.
Entrance constraints are not working/prefix bounds.

## Physical gauge changes are not free uncertainty deletion

Rebasing position by a constant p_0 while retaining S_dot=p entails

    p' = p-p_0,
    S' = S-S_0-(t-t_0)*p_0.

A fixed zero S pseudo-observation is not invariant under this transformation
unless its target/forcing is transformed as well. Therefore simultaneous
claims of p_true=S_true=0 at handoff cannot silently retain the old bounded
physical primitive and the old forcing ledger. An S-only integration-origin
choice likewise requires an explicitly stated source/reference convention.

## Remaining mathematical obligations

This entry lemma does not create a COMPLETE-BRMM source-uniform cover or a
certified maximal basin. Numeric physical velocity/position/primitive bounds,
or their declared correlated forcing alternatives, must be qualified. The
same-source frontend/tuner/scheduler and reachable P/H/R relations must be
propagated jointly before constructing compatible word storage, coercivity,
endpoint contraction, every-literal-prefix bounds, outward augmented LDLT
and first-exit retention. Finite-precision charges must be consumed there.

BIAS0, BIAS1 and BIAS2 retain their own physical-driver recurrences and the
family-parametric projection/Joseph graph. Positive BIAS2 separation is not
introduced as a prerequisite merely for bounded bias plus practical motion
ISS. Full startup capture follows only after genuine P4 closure, with the
real timeout, magnetometer, tuning, covariance and reset branches retained.

Current issues are entry-model/dependency/source/capture obligations (D/C/E/F),
not evidence of an admissible nonlinear instability (A). No uniform capture
time, certified maximal basin, endpoint contraction or indefinite stability
is claimed by this document or its exact entry-identity tests.
