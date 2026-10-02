# OU-III stability proof handoff — PR #641

## Controlling path

construction/capture -> magnetically informed H18 -> H18/LIN BIBO ->
compact H18/A21 release -> retained outer A21 -> finite outer-to-inner entry ->
V<0.15^2 local LaSalle -> linked finite-error/prefix retention ->
nonlinear/transition/float32 closure.

Do not revive the retired scalar kernel-ceiling, independent dH, arbitrary
covariance-box, AW-only passivity, pointwise global AW tracking, or
zero-homogeneous-action => zero-base-innovation routes.

## Physical qualification used in the current conditional branch

MARINE candidate:
- T_E=T_P=30 s
- theta_E=2 degrees
- P_E=0.03 m

SLOW+FAST candidate:
- H_a=60 s, C_a=0.05 m/s
- H_g=60 s, C_g=0.002 rad

These FAST values are candidate assembled-device qualifications, not measured
specifications. The theorem remains conditional on their independent device
qualification.

## Closed / retained

- Local V<=0.15^2 theorem retained. P_aw,aw<=16.48 I and
  |a_hat_w-a_phys|<4.06 sqrt(V); 17-s field-axis strict margin
  0.10527117647 m/s^2.
- Historical fast-gyro periodic witness needs about 0.04 rad half-cycle signed
  accumulation and is excluded by candidate C_g=0.002.
- 17-s continuous SLOW+FAST gyro field-axis charge:
  0.0136432352941176 m/s^2.
- Real-arithmetic source-only A* charge conservatively:
  1.0794348403 m/s^2 < g/5=1.96133, margin 0.8818951597 m/s^2 before float32.
- Quaternion small-angle polynomial is already part of qualified gyro
  prediction transport; finite d^2 v reset curvature belongs to linked
  finite-error supply, not source-only AW forcing.
- Complete regular-word joint zero-action kernel has nullity <=1:
  four-S/process removes independent LIN/AW root; MAGNETIC SERVICE reduces
  root attitude image <=1; active A21 J_ba=I kills pure BA and makes BA unique
  for a surviving attitude amplitude. Additional rows may kill the line.
- Persistence of a surviving physical tilt/BA compatibility line through
  qualified MOVING superwords is excluded by A*. Compactness therefore gives
  finite-superword homogeneous strict dissipation at existence level.
- Captured H18 release is compact. Recurring root covariance lower bound makes
  V bounded on the release image, hence a finite compact outer A21 region
  exists. Compactness also gives a finite practical absorbing radius.

## Important non-claims

- The old explicit G0 numerical floor is NOT promoted: its stronger
  m_perp<=0.4 and nominal L1 force<=1.2 premises are not source-uniformly
  proved by the new A* bridge.
- No numerical finite-superword contraction constant is claimed from
  compactness alone.
- Outer-to-inner entry below 0.15 is NOT yet proved.
- End-to-end theorem, prefix retention, float32 and regime composition are
  NOT closed.
- The retrospective linked defect e_N-M e_0 is diagnostic only and is marked
  ineligible for theorem entry.

## Immediate blocker

For e_N=M e_0+b use the exact completed-square quantity

    G_gamma = J0 - M' JN M - gamma J0
    z       = M' JN b
    chi     = b' JN b + z' G_gamma^-1 z

The target is

    sup_h inf_gamma chi_gamma(h)/gamma < 0.0225

on the SAME reachable history, plus corresponding every-prefix inequalities.

The true forcing b_W must be composed from literal local defects:

    d_k = e_{k+1} - A_k e_k
    b_{k+1} = A_k b_k + d_k.

Proof infrastructure now exports physical p,v,S,a, quaternion, BG, BA and
estimator xext/qref at internal literal boundaries, reconstructs the 21-state
error, and hard-gates the carried observer on

    ||b_local - (e_N-M e_0)||_inf < 5e-10.

The next conversation should first make that native gate green for quiet and
wave, then feed b_local into the optimized linked ratio. Only after carried
feasibility passes should it build a source-uniform SLOW+FAST/prefix enclosure.

## CI / evidence hygiene

Proof-driver instrumentation changes invalidate committed diagnostic
fingerprints. Regenerate affected evidence; do not weaken fingerprint checks.
Do not merge theorem flags beyond the mathematical statements actually closed.
