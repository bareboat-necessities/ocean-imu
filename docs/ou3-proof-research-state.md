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
