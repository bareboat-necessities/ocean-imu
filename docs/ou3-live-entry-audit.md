# OU-III executable Live-entry and shared-source error graph

## Scope and status

The source of truth is the shipping `SeaStateFusion_OU_III` outer wrapper,
`SeaStateFusionFilter_OU_III` inner wrapper and `Kalman3D_Wave_OU_III` core.
`tests/kalman_ou_iii/live_entry_audit-test.cpp` executes them without replacing
filter code. `ou3_p4_live_entry_graph.py` checks the exact shared-origin algebra
and carries the continuous physical jet and signed source-energy transform.

P4 motion and P5 capture are **not closed**. BIAS0, BIAS1 and BIAS2 retain their
separate closed conditional driver/projection lemmas. P3 delta is `1e-18`.
No shipping header, physical source, quality threshold or deployed arithmetic
is changed by this audit. The old independent entry product is retained only
for reproducible conditional diagnostics, not treated as startup reachability.

## The eight initialization questions

1. **Is estimated velocity zero at fresh Live entry? Yes.** The core constructor
   zeroes it, and the normal outer pre-Live path does not run core motion
   prediction or measurement correction. `goLive()` does **not** zero it.
2. **Is estimated position zero? Yes, for the same fresh path.** There is no
   additional position reset at the Live transition.
3. **Is estimated S zero? Yes, for the same fresh path.** It is carried from
   construction; it was not integrating the warmup position uncertainty.
4. **Where are the actual resets?** Fresh `begin`/inner `initialize[_ext]`
   constructs the zero core. `initialize_from_truth` is a separate test/simulation
   API that loads p/v/aw and zeros S and biases; normal startup does not call it.
   Attitude initialization and tilt re-lock reset the attitude error coordinates,
   not the additive navigation or bias values. H18/A21 release is not a reset.
5. **Once or propagated?** Fresh-path zeros persist until the first core event
   after `goLive`. Directly driving the inner filter earlier produces a different
   reachable state; the native sentinel test confirms that nonzero v/p/S/aw and
   biases survive `goLive`. No universal zero claim applies to that API history.
6. **Covariance?** The fresh v/p/S diagonal blocks are I, 400I, 2500I, with zero
   mutual cross-covariances. Attitude initialization rebuilds its covariance,
   clears attitude/bias cross terms and the entire BASE_N-by-12 base/linear
   block (including gyro-bias/linear terms). Live replaces the aw covariance
   with the same active stationary OU covariance and clears all aw cross terms.
   H18-to-A21 preserves mean states and motion P, flooring BA diagonal entries.
   None of these covariance choices establishes a hard true-error bound.
7. **Truth origin?** No executable handoff resets physical truth. The paper's
   session reference is S_true(t)=integral_0^t p_true(s)ds, with p describing
   physical NED displacement. A nonzero wave displacement at handoff is not
   renamed zero by setting a filter estimate to zero.
8. **True errors?** At the fresh boundary, before subsequent completed events,
   e_v=v_true(t_L), e_p=p_true(t_L), e_S=S_true(t_L), e_aw=a_true(t_L),
   e_bg=bg_true(t_L), e_ba=beta_true(t_L). These share one physical/startup history.
   Neither zero truth velocity nor zero truth position is established.

## Other state transitions

| Event | State values | Covariance | Frontend/scheduler |
| --- | --- | --- | --- |
| Fresh construction/Cold | q initially identity; additive means zero | Constructor blocks; BA held/decoupled | Frontend reset; pseudo elapsed=0 |
| Proxy learning/TunerWarm/Ready | Proxy and adaptation evolve; core additive means held | No core Riccati prediction on outer pre-Live path | WPE, band, noise covariance, moments, candidate/active commits evolve; no pseudo tick |
| First `goLive`/H18 | Proxy attitude seeds core; all additive means preserved; Mahony integral is **not** copied to core bg | Attitude/cross reset; aw stationary reseed | Learned tuner history committed, not restarted; period retarget preserves elapsed progress |
| First subsequent Live sample | Actual prediction and due S/accelerometer events | Full reachable Riccati/Joseph/reset updates | Actual held deadline and candidate/active ordering |
| H18 to A21 | Mean values preserved | BA enable floors diagonal; motion P retained | Magnetic/hold guards decide release |
| Live tilt re-lock | Yaw-preserving attitude reinitialization; additive values preserved | Attitude/base-linear cross terms cleared as in initializer | Hybrid guard/cooldown retained; not a new zero-motion entrance |

The core has a full 21-state covariance while BA is held; H18 denotes the
motion-error subsystem, not allocation of a new 18-state filter. Repeated
initialization, configuration overrides, callbacks after the handoff sample and
hybrid resets must be represented separately, not hidden inside the fresh graph.

## Why the independent 300 m*s factor is not the actual entrance

The retained history shows 300 as a declared handoff budget in the operating
contract (already present at its `a42cc67` relocation into `tools/stability`)
and copied into the closure domain by `c8e54b7`. Neither introduces a derivation
from the shipping startup. It is not justified by the constructor's 50 m*s
covariance standard deviation. Earlier provenance must be followed across the
file relocation; the reproducible audit workflow retains that history.

The native test uses a C2 periodic physical continuation: every p coordinate
is +2.1 m for 200 s, transitions by a quintic smoothstep to -2.1 m, holds,
and returns, with an 800 s cycle and half-cycle antisymmetry. All derivatives
match at joins. Its velocity, acceleration and integral primitives are globally
bounded, it has no DC displacement or acceleration over a full cycle, and its
vertical Hs is at most 8.4 m. The exact rational certificate bounds every
continuous 10 s window, not only sampled windows. With the audited candidate
quiet threshold 0.009 and impulse cap 2, every such window is Q. This is an
analytic witness within the parameterized bounded-primitive BRMM declaration
when its displayed constants are admitted, not a replacement source family.

On the initial plateau the true acceleration is zero. The unchanged wrapper
with quiet IMU input reaches Live at sample 18051 (physical time 90.255 s):
S_hat=0 but the session-origin S-error norm is about **328.285 m*s**.
With the bounded zero-proxy sensor sequence described below, the same physical
plateau gives **545.632 m*s** at the timeout handoff. A smaller arbitrary S
radius would not repair this false entrance assumption.

The old frozen-word Cayley value 4.5788 exceeds its own chart-valid enclosure
range. It is not a rigorous nonlinear chart-exit counterexample. The values
0.9411/0.6537 obtained by omitting the independent factor are also diagnostics,
not endpoint/prefix certificates for the corrected source graph.

## Exact shared-origin lemma: a wider correlated entrance

At fresh entry, e_S=S_true, so the actual regularizer residual is

    r_S = -S_hat = e_S - S_true = 0.

Set c=S_true(t_L) ONCE, for the entire subsequent Live history, and define

    S_L(t)=S_true(t)-c,     e_S,L(t)=e_S(t)-c.

Now S_L(t_L)=e_S,L(t_L)=0, while r_S=e_S,L-S_L remains exactly unchanged.
This is a proof-reference transformation, not a new filter, a truth reset in
shipping code, or an altered BRMM source. For arbitrary same-event K,

    (I-K H_S) E_S + K = E_S.

Prediction preserves E_S; non-S observations annihilate it. Thus a common
physical-S/error-S origin column remains exactly E_S through **every completed
literal event prefix** and is removed by the centered error row operation.
All finite residuals, quaternion corrections, reset transports, projection
branches, P/H/R/K, tuner state and deadlines are identical. Because the origin
is never an executable input, binary32 execution is bit-for-bit identical too.
The producer checks the formal matrix-polynomial identity in H18 and A21.

There is consequently no attitude-chart limit on this shared constant-origin
direction. This is stronger and more physical than shrinking an independent
S ball. It does NOT delete the three evolving S-error coordinates from the
18-motion metric or re-zero S at later word boundaries.

Working/prefix regulation remains necessary. The same continuous physical history
gives, for h after handoff,

    S_L(h)=h*p_L+h^2*v_L/2+integral_0^h (h-s)^2*a(t_L+s)/2 ds.

The acceleration moments, p_L and v_L are correlated source coordinates. A
consequent norm bound is min(h*P_m, h*|p_L|+h^2*|v_L|/2+A_m*h^3/6, 2*S_m),
using only terms whose uniform primitive bounds have actually been admitted.
Likewise the source-energy transform must keep its signed cross term:

    integral |S-c|^2 = integral |S|^2 - 2*c.integral(S) + T*|c|^2.

A certificate cannot remove the entry factor and silently drop this supply.
Recover original-reference errors with |e_S|<=|e_S,L|+|c|. An original-reference
uniform ultimate bound still needs a uniform bound on c. Coercivity of a
centered metric is not automatically coercivity about the uncentered origin.
Position cannot be removed in the same way: shifting p by p_L adds the secular
term -h*p_L to S and can violate the indefinite bounded-primitive premise.

## Actual Live can retain the prior forever on quiet input

The quality path requires tuner readiness; the aligned timeout path does not.
In particular it does not require a usable WPE period or north_ready. The
native test constructs a bounded scalar accelerometer disturbance, at most
0.03326797485 m/s^2 on the tested host, that cancels the private Mahony proxy's
normalization residual on a level, zero-rate physical history. It executes the
shipping observer and verifies its up signal is exactly zero at every startup
sample. Both magnetic and nonmagnetic test configurations enter runtime Live
at sample 30002, shipping clock 150.0007477 s, with WPE usable=false.

For a zero-motion continuation, zero input and zero WPE filter/moment states
are an invariant set: no positive variance ratio exists to establish a period.
The bounded sensor witness therefore invalidates a universal finite-capture
claim into a measured-period-only domain when that disturbance is admitted.
It is not filter instability, and expiry of the timeout is not a certified
T_capture. Host binary32 reproduction is not target-toolchain qualification.
With exactly nominal quiet sensor samples this host instead creates a small
proxy transient and enters with a usable period at 90.245 shipping seconds;
that empirical behavior must not be assumed for every disturbance realization.

The proof-side joint transition now retains the actual 0.2 Hz prior, the
uninitialized log state, first-valid initialization, and the one-way measured
period takeover. The current tuner uses the PREVIOUS usable selector; its band,
variance/noise state, candidate/active EMA and commit remain in the same history.
No independent frequency/sigma rectangle is introduced.

Two existing interval guards also silently dropped possible valid updates when
variance or omega^2 straddled its threshold. Both sides are now emitted. The
common moment weight cancels exactly in the ratio: var_v/var_eta=C_v/C_eta.
Only strictly qualified post-moment-start cells are accepted; an unhandled
weight guard raises, rather than certifying a subset. Binary32 moment-cancellation
error remains a separate required enclosure, not covered by the real identity.

## Remaining certificate and limiting boundaries

No maximum numerical P4 basin has been certified. The exact shared S-origin
direction has no correction/chart restriction, but the other boundaries still
require a joint reachable cover: attitude/gyro bias from the actual proxy;
v/p/aw from the same physical history; BA from its own family and radial
projection; and tuner/P/scheduler from the actual startup and transitions.
The numerical BRMM V_m, P_m, S_m and recurrence/excitation budgets remain
unfrozen in `ou3_brmm_contract.py`. Do not infer their values from Hs or covariance.

The missing source-uniform cover must include prior/measured adaptation,
Q/O/mixed motion, all BIAS histories, physical S forcing, reachable P/H/R and
same K, every reset/projection and literal prefix. Source-dependent compatible
storage, coercivity, rho_H<1, Gamma, C/C_p, first-exit retention and all domains
must then be checked by outward augmented LDLT. The current producers do not
supply those matrices. No local lemma or point diagnostic closes them.

Failure classification: **D** for the independent fresh-entry model; **C** for
lost WPE branches and enclosure primitive boundary failures; **E** for missing
source qualification/full cover; **F** for capture into a measured-only subset
on admitted quiet disturbances. No admissible nonlinear instability (**A**) is
established. A claim of exhaustive infeasibility of all compatible-storage
constructions (**B**) would also be unsupported.
