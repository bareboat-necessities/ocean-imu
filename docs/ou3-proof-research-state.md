# OU-III proof: current continuation state

Source audit: PR #637 at `1bbd3b3bc58879d7f95863cef5726ee10c80616f`.
The complete preceding ledger is preserved byte-for-byte in
[the pre-field-axis-audit ledger](ou3-proof-research-state-before-field-axis-audit.md).
The long-form derivations remain in `ou3-corrected-word-proof.md` and
`ou3-linked-soft-return.md`; do not revive their explicitly retracted claims.

## Single proof path and unchanged scope

construction -> capture -> magnetically informed H18 -> refinement/release ->
recurring A21 -> regional practical stability. The contraction word remains
100 s; 17 s remains nuisance/root warm-up. Runtime, physical assumptions,
coupled tau/sigma_aw/r_S/T_S chronology, noise, and quality gates are unchanged.
Stillness is not removed from the physical class by this continuation.

## New analytical result

[Source-audited regional field-axis exclusion](ou3-field-axis-regional-exclusion.md)
(FA1--FA12) is the current calculation. It repairs two source mismatches in QR:
accelerometer updates occur every regular IMU sample (maximum qualified applied
gap 6 ms, not 40 ms); default AW synchronization is an additive PSD floor, not
congruence. Isotropy makes its marginal eigenvalues exactly max(old,target),
which supplies the needed ceiling without a full 21-state covariance ceiling.
An exact-rational audit of the literal Q polynomials supports the conservative
regular real-arithmetic bound P_aw,aw<=16.48 I, hence ||e_aw||<4.06 sqrt(V_base).

The sharp same-history trapezoidal inequality gives

    sum_i w_i ||P_bi(a_hat_i-g_model)||^2 >= max(M_ref(r),0)^2,
    M_ref(r)=G0-delta_g-eta_ref(Amax+g_model_max)
             -4.06(1+eta_ref)r-2Vmax/L-Jmax h_acc/4.

Every a_hat_i and V_base,i is evaluated at the actual PRE-ACCELEROMETER prefix.
For a fixed, correctly bounded nominal reference with inclination <=80 deg,
the 100-s pre-radius margin is 1.4429069015... m/s^2. At r=.25 the margin exceeds
.4279 m/s^2 before reference-error charges (.4278 after worst-case word-boundary trimming). This is a conditional local
exclusion of the particular field-aligned base trajectory, not stability.

## Scope correction / unresolved implication

The physical magnetic field and the learned nominal reference are different
objects. MagAutoTuner's initial 5% horizontal-fraction gate gives a weaker
positive fixed-reference margin, but the default continuous hard-iron path
rewrites the reference without reapplying that fraction gate. This does NOT
prove it fails in a shipping execution. It prevents treating the startup cone
or physical 80-degree cone as an already established all-time nominal cone.

Next: prove eta_ref (or a sufficient direct nominal cone/variation estimate)
from literal acquisition, refinement and hard-iron updates, then discharge
prefix retention in the same physical history. Do not invent zero defects,
freeze the default reference, disable continuous hard iron, replace actual
accepted-update times, or assume stability while proving retention.

## Classification and validation

PROVED CONDITIONAL: FA1--FA12 under their stated regular covariance, prefix,
reference and physical premises. The new bound is a base geometric row-energy
floor, not the full innovation-weighted nuisance-eliminated word action.

OPEN: unconditional shipping exclusion, capture/release, source-uniform linked
O2 return bound, numerical full-word rho<1, nonlinear/physical supply retention,
and target arithmetic. No theorem or quality flag is promoted. Zero homogeneous
action must never be substituted for zero BASE innovation. Regional r is a
full-state covariance-metric radius, not the retained six-degree attitude angle.

The 13 new exact scalar/rational regressions pass with
`python -m unittest discover -s tests/validation -p 'test_ou3_field_alignment_exclusion.py' -v`.
These validate the analytical substitutions, not a full trajectory enclosure.
Full native validation and CI were not run here. Inherited failures, numerical
experiments, unsuccessful routes, and exact previous status are retained in the
archived ledger and previous PR history; they are not claimed repaired.


## Shipping-reference and prefix-premise audit — 1c8f4e94

The requested follow-up is now source-audited in
`docs/ou3-shipping-reference-prefix-retention.md`.

The continuous hard-iron reference update has an exact same-statistics
Lipschitz law: with `L(b)=wbar-Abar*b` and `||Abar||<=1`, each canonical
horizontal/z reference component moves by at most the applied body-bias
increment, and the full canonical reference by at most `sqrt(2)` times it.
The 45-s slew gives the corresponding per-update increment. This is literal
shipping chronology, not a compatibility relaxation.

That result does NOT yield the needed all-time nominal cone from the present
contract. Startup/refinement enforce a 5% horizontal fraction, but the default
continuous path does not reapply that fraction gate; it only requires positive
horizontal magnitude above .001 uT. Its loose accepted-fit envelope
`.35*(75+5+2)=28.7 uT` exceeds the declared physical 15-uT horizontal
minimum. More decisively, its statistics use the private Mahony tilt, and no
all-time deterministic true-to-proxy tilt tube has yet been proved. Therefore
the physical-field cone cannot be silently transferred to the nominal
reference.

The FA12 full-storage premise also cannot be inherited from six-degree capture.
The fixed-reference field-exclusion radius must satisfy `r<.355396`, whereas
the handoff tilt covariance sigma is .035 rad, so a six-degree tilt error alone
has minimum covariance-metric radius `(pi/30)/.035=2.99199...`. This does not
say handoff has six-degree error; it proves the existing capture target does
not imply the small FA12 storage ball.

To remove that artificial circularity, the field-axis lemma is reformulated in
the exact component it needs. If
`||a_hat_w-a_phys||<=eps_aw` at every relevant prefix, the fixed-reference
100-s exclusion needs only
`eps_aw < g*sin(10deg)-.26 = 1.4429069015... m/s^2`.
The varying-reference version is given as PR10 in the new note. This is much
weaker than requiring the entire 21-state Mahalanobis error to be <=.25.

Prefix invariance remains a simultaneous fixed-point problem, not an upstream
premise: the same block must close kernel covariance return, word error supply,
and every-prefix supply. Proving `V<=r^2` first from covariance bounds would
be circular because the strict word return needed for that storage recursion is
the open O2/BP obligation.

Next decisive calculation: derive (or falsify) an all-time private-Mahony tilt
tube and an all-prefix AW tracking tube from their literal error equations on
the SAME physical history, retaining the coupled tau/sigma_aw/R_S/T_S
chronology; insert those component tubes into BP-4/BP-10 and solve the retained
rectangle simultaneously. No theorem flag is promoted.


## Correction: the 1.4429069-m/s^2 target is sufficient, not yet implied — current head

The immediate blocker is still ONLY the persistent nominal force/field
collinearity trajectory. The desired contradiction is

    |weighted mean P_B(a_hat_w-a_phys)| < 1.4429069015... m/s^2

for the 100-s, 6-ms, 80-degree specialization, because physical bounded
velocity+jerk contributes at most .26 m/s^2 while the transverse gravity
requirement is g sin(10 deg)=1.7029069015....

However, the existing sections 59--61 of
`docs/ou3-aw-adjoint-cancellation.md` already prove that this inequality
CANNOT be inferred merely by saying that a_hat_w, tau, sigma_aw, R_S and T_S
are coupled.

For a prescribed smooth periodic physical history, the measurement-only
front end and tuner determine one periodic applied coefficient word U_*.
Conditional on U_*, covariance propagation is independent of innovation
VALUES. After the stabilizing periodic Riccati orbit is fixed, the mean is an
affine periodic linear system. The finite Kalman map from a measurement word
to its innovation word is block lower triangular with identity diagonal and is
therefore invertible. The exact periodic compatibility equation is PR9/PR14.
Generically nonsingularity gives a UNIQUE compatible forced periodic orbit; it
does not exclude one.

Therefore the coupled tuning law is essential for fixing the coefficients, but
it is not an amplitude theorem for deterministic estimation error. Covariance
P_aw,aw<=16.48 I likewise bounds uncertainty/action geometry, not
|a_hat_w-a_phys| on an arbitrary deterministic forced execution.

This invalidates the proposed shortcut

    coupled tuner law => |mean(a_hat_w-a_phys)|<1.4429.

The 1.4429 number remains a correct sufficient threshold. To prove the
pathological trajectory inadmissible one must instead show that the UNIQUE
self-consistent PR14 solution violates an EXISTING admissibility condition
(physical p/v/a/jerk, bias/projection, gate, service, retained local angle), or
prove a sharper signed identity that forces such a violation.

The constructive audit has already made this falsifiable. The first
commensurate 6-s candidate is excluded because its literal late relative
attitude reaches 7.488 degrees >6 degrees. The analytically refined 12-s
candidate has pre-compensation physical amplitudes strictly inside the
declared envelopes and is the decisive next target. It has not yet been
certified as a shipping counterexample because its late periodic
tuner/covariance/mean orbit, six-degree bound, gates/service and complete
weighted functional remain to be enclosed.

Accordingly, do NOT claim the accel||mag pathology is excluded yet, and do NOT
spend the next calculation deriving a generic AW tracking tube from covariance.
The decisive calculation is the literal 12-s PR14 periodic orbit with outward
enclosure. If it violates an existing condition, extract that violation as the
analytical exclusion lemma. If it satisfies all conditions, the pathological
trajectory is admissible under the current theorem contract and the stability
proof must be reformulated; no additional physical assumption may be silently
introduced.


## Controlling LaSalle update — radius-local field-axis candidate excluded

The controlling local proof is now `docs/ou3-radius-local-field-alignment.md`.
It replaces the attempted global nominal-AW/innovation control for the
zero-dissipation invariant-set question.

Inside the retained storage ball, the current literal source audit gives
`||a_hat_w-a_phys||<4.06 r` from `P_aw,aw<=16.48 I`; the inherited BA
marginal gives `||e_ba||<=r/40`. BA is NOT added to the field-alignment
tube because the implemented accelerometer attitude Jacobian uses the nominal
CoG vector `a_hat_w-g`; BA is a separate measurement column and lever arm is
attitude-independent there.

On the explicit committed-field branch `sigma_w>=1/5`,
`||P_B g||>=9.80665/5=1.96133 m/s2`. For one MARINE history,
`||v||<=5.5` gives a continuous T-window point with
`||P_B(a-g)||>=1.96133-11/T`. The 100-m/s3 jerk bound and an ACTUALLY
APPLIED accelerometer gap <=.006 s transfer this to a sampled epoch with

    m_phys(T)=1.96133-11/T-.6.

At T=17 s, `m_phys=0.714271176470588...`. Hence persistent nominal
field alignment is impossible whenever
`4.06 r<m_phys`, i.e. `r<0.175928...`. Adopt the deliberately
conservative local exclusion radius

    r_FA=.15, T_FA=17 s,

with strict margin `0.105271176470588... m/s2`.

Conditional on the already established zero-dissipation classification
(leaving only the field-axis candidate), this proves

    Inv_MARINE({D=0}) intersect {V<=.15^2} = {0}

on the regular real-arithmetic retained branch. The existing closed-stratum
compactness argument then gives existence of finite `m` and `eta_D>0`
for homogeneous finite-window strict dissipation by contradiction/diagonal
extraction. No numerical eta_D is claimed.

OPEN after this local invariant-set closure: entry/every-prefix retention in
V<=.15^2, finite capture/H18/release, nonlinear/source supply over the
finite-window block, recurring transition budget, and full float32 totality.
The end-to-end regional practical-stability theorem is NOT promoted.

Historical O1/O2/kernel-ceiling and signed-reader calculations remain useful
research but are non-controlling for exclusion of the zero-dissipation
field-axis trajectory.


## Entry/retention continuation — staged A21 entry required

`docs/ou3-staged-entry-retention.md` is now controlling downstream of the
r_FA=.15 field-axis exclusion.

Direct H18 release into full V<=.15^2 is NOT a valid universal target.
Shipping held BA is decoupled with sigma_bacc0=.004 m/s2; enabling A21 merely
floors its diagonal variance to sigma_bacc0^2. The admitted physical BA norm
is .22516660498395405. A permitted held b_hat_a=0 therefore has decoupled
release contribution V_ba=(.22516660498395405/.004)^2>3168, versus .0225 for
the final local ball. No assumption or runtime change is made.

Homogeneous every-prefix retention IS closed once inside the local ball:
prediction/correction are covariance-metric nonexpansive and congruent resets
preserve storage; PSD covariance inflations cannot increase fixed-error
storage. Finite nonlinear/source/arithmetic residuals require an inner root
radius r_in<.15 with

    sup_k G_k(.15,d) <= .15-r_in,
    E_W(.15,d) <= (1-q) r_in,

q=sqrt(1-eta_D).

The correct entry path eliminates the held BA coordinate. At release it is
literally decoupled; after release use the BA Schur complement / quotient
storage V_o|ba. Prove quotient-local entry/retention at .15, apply the
field-axis LaSalle exclusion there, then use active A21 BA dynamics and the
dissipative projection sector to enter the full ball. This is the same
H18-complement architecture, not a new proof path.

Current limiter: source-uniform outer/BA-quotient release set inclusion and
quantitative quotient finite-window dissipation/supply. Do not retry universal
direct full-V release.


## Literal release outer-storage audit — attitude is the entry obstruction

A dedicated unchanged-header release snapshot now evaluates the carried
construction history at first A21 activation (step 36008, 180.039996 s).
Correct BA elimination uses the outer covariance marginal; at the release
boundary BA cross covariance is zero so it also equals the conditional outer
block there.

Release tilt is 8.144927 deg. The attitude covariance eigenvalues are
3.5352043e-6, 3.5550410e-6 and 1.2944541e-5. Minimizing the BA-eliminated
outer storage over every other non-attitude outer coordinate still gives

    V_outer,elim >= theta' P_theta^-1 theta
                  = 5716.295063726523,
    sqrt(V_outer,elim) >= 75.60618403098071.

Thus this literal release is nowhere near r_FA=.15; ATTITUDE alone consumes
the radius. At this covariance, V<=.15^2 would necessarily require tilt
<=.000539678 rad=.0309213 deg.

This is finite carried evidence only: the stress history does not certify
all-time MAGNETIC SERVICE and therefore is not an eventual-capture
counterexample. It does prove that stage flags/release mechanics themselves
do not imply entry into the tiny local storage ball. The controlling proof now
needs an outer A21 retained/capture region and finite entrance from that region
to the inner r_FA=.15 LaSalle ball. The r_FA field-axis inequality cannot
simply be enlarged to the observed release storage: its current source-audited
AW component conversion loses positivity above r~=.176.


## Literal release outer-entry calculation — AW is the blocker

The requested BA-eliminated release calculation is recorded in
`docs/ou3-release-outer-entry-audit.md`.

A dedicated unchanged-header diagnostic snapshots first A21 BA activation at
step 24016 (120.079997316 s) on the carried diagonal-wave construction. Using
the exact physical construction and outer marginal identity
`V_o=min_ba V=e_o'P_oo^-1 e_o`, the finite carried release has
`V_o=79300.039948`, sqrt=281.603, versus target .15.

Principal-block lower bounds after minimizing every other coordinate identify
AW as dominant: attitude 9.7976, bg 3.59e-7, v 331.982, p 1996.290,
S 84.977, AW 60867.918. The release AW error norm is about 8.385 m/s2.
This finite history is not an all-time service/capture certificate.

The second requested calculation has a negative but decisive result:
the present assumptions do NOT imply source-uniform pointwise AW entry.
An existing admitted A21 history already refutes pointwise physical-AW
tracking at 7.647 m/s2 on a 16-s window, whereas the r_FA=.15 local lemma
requires <.609 m/s2. P_aw,aw<=16.48 I only converts an ALREADY SMALL storage
to an AW component bound; it does not bound deterministic AW error before
entry. No AW mean projection supplies such a cap.

Therefore direct H18/release -> V_o<=.15^2 is not the correct bootstrap.
The local FA theorem remains valid after entry, but entry must use a shaped
signed/windowed AW functional or return to the global same-history FA
reachability/action calculation. Do not retry a source-uniform pointwise AW
tracking lemma; it is already falsified on the admitted class.
