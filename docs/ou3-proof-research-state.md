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
