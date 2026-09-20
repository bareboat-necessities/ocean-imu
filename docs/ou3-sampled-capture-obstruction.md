# The stated motion domain does not imply six-degree physical capture

The controlling failure is now a proved obstruction to the **capture and
retention** part of the requested theorem. It is not a failure of the nominal
Kalman energy identity, a claim of divergence, or a quiet-water exclusion.
Both histories below have nonzero motion. They satisfy the declared all-time
motion, bias and point-sampled sensor bounds. The common regular
real-arithmetic shipping execution also satisfies the actual one-second
magnetic-service inequality, with its carried covariance and corrections.
Nevertheless its physical tilt error stays at `acos(3/5)`, about 53.13 degrees.

The six-degree target is the current
`a21_retained_domain.tilt_error_max_rad` in `constants.json`. No physical
assumption, shipping operation, or quality threshold is changed here.

## 1. Two continuous physical histories, with identical sensor samples

Use body-to-world attitude and put `g=9.80665`, `h=1/200`, `w=2*pi/h`,
`c=3/5`, `s=4/5`, and world magnetic field `B=(75,0,0)` microtesla. For
`epsilon=+1` or `-1`, define one persistent history by

```
R_e = [[1, 0, 0], [0, c, -e*s], [0, e*s, c]],
A_e = g*(0,e*s,1-c),
a_e(t) = A_e*cos(w*t),
v_e(t) = A_e*sin(w*t)/w,
p_e(t) = -A_e*cos(w*t)/w^2,
q_e(t) = -A_e*sin(w*t)/w^3.
```

Exactly, `p'=v`, `v'=a`, `q'=p`; R is constant and body angular rate is zero.
The physical gyro and accelerometer residual biases and sensor residuals are
zero. In particular, no estimator OU prior is imposed on physical truth.
Neither q nor any other physical coordinate is restarted at a proof boundary.

The all-time bounds are

`||a|| <= 2g/sqrt(5) < 8.8`,
`||v|| <= 2g/(sqrt(5)*w) < .00732 < 5.5`,
`||p|| <= 2g/(sqrt(5)*w^2) < .00000610 < 8.1`,
`||q(t2)-q(t1)|| <= 4g/(sqrt(5)*w^3) < 1.02e-8 < 1100`.

The rational certificate uses only `pi>3` for these conservative bounds.
The displacement has zero DC; its primitive is globally bounded. These are
not independently chosen acceleration samples with incompatible kinematics.

At every IMU time `t_k=k*h`, `a_e(t_k)=A_e`, and exactly

`R_e^T (A_e-g*e_z) = -g*e_z`,
`R_e^T B = B`, `gyro=0`.

Thus **both nonzero histories produce the same entire sample sequence**:
accelerometer `(0,0,-g)`, gyro zero, and magnetometer `(75,0,0)`. Magnetic
callbacks every eight IMU samples use the deployed 25 Hz cadence. Their norm
and horizontal component satisfy the field limits. All samples are finite,
unsaturated and have zero physical measurement residual in the stated
point-sample model. Temperature can remain fixed at its compensation reference.

A proven sensor integration/anti-alias model could exclude this construction.
The present formal premises do not contain one. This is a counterexample to
that formal domain, not a claim that a vessel normally executes a 200 Hz wave
or that an unmodelled real sensor has a flat frequency response at 200 Hz.

## 2. Literal shipping construction, startup and release

The common input initializes the proxy level. Its gyro and gravity correction
remain zero, so its attitude is I. The world gravity average is `(0,0,-g)`;
after its finite warmup and hold the aligned-branch and gravity-quality
predicates hold. Constant valid north samples pass the norm and horizontal
MagAutoTuner gates and complete the finite acquisition count/window. The
learned reference is `(75,0,0)` and the yaw gauge is zero.

Zero wave input need not qualify the measured-period gate. This does not
block this history: the literal handoff timeout requires proxy readiness and
the aligned gravity branch, which already hold. It therefore hands off no
later than the configured finite timeout, with a gauged level attitude.

After handoff the nominal attitude remains I, and the gyro bias, LIN means
and accelerometer bias remain zero. Each prediction preserves these values;
S, acc and mag innovations are zero. Attitude injection/reset is identity and
bias projection is inactive. Parameter adaptation can change covariance,
noise weights and S cadence, but multiplies zero nominal LIN means.
The reference-refinement window also receives constant valid north in a
level proxy frame. It completes in finite time, sets the same field/yaw,
and releases the external BA hold. The actually applied magnetic count and
one-second guard complete the internal release. Constant proxy tilt has zero
excitation span, so the continuous hard-iron fit cannot qualify a new offset. Hence A21 is reached in
finite time on this very execution, without an artificial covariance reset
or a seeded post-startup replacement trajectory.

The native regression `sampled_capture-test` exercises these predicates on
the literal float shipping wrapper. It verifies finite Live/refinement/A21
entry, every subsequent supplied magnetic correction being applied, a level
attitude, and exact separation of the heading covariance block. This finite
regression supports the source binding; the all-time argument is the
invariance above and the covariance/service proof below, not extrapolation
of that replay.

## 3. All-time actual magnetic service, not packet cadence alone

At this stationary nominal state, `(theta_z,b_gz)` is an invariant covariance
and linear differential block. S and accelerometer sensitivities have no
column in this pair. The two nonzero magnetic rows separately observe
`theta_z` and `theta_y`; all cross-axis covariance is initially zero and is
preserved by the source prediction, Joseph corrections, AW sync and identity
attitude reset. This is an exact symmetry reduction, not an information lift
from an arbitrary principal block.

After dividing the magnetic heading row by its field magnitude, the pair is

`F(t)=[[1,t],[0,1]]`, `H=[1,0]`, `r=(.8/75)^2`,

`Q(t)=[[sg^2*t+qb*t^3/3, qb*t^2/2], [qb*t^2/2, qb*t]]`,

where `sg=.00135`, `qb=1e-10`. At zero rate the literal source rotation and
Simpson integral give exactly these matrices. Q is positive definite; the
source 6x6 LDL hygiene does not add a floor in real arithmetic.
The eight IMU predictions between magnetic updates compose exactly to
`t=.04`. Use the two covariance bounds

`U=diag(.008,.00002)`, `U_prefix=diag(.009,.00003)`.

The gauged handoff gives `P_theta=.087^2<.008` and `P_bg=1e-6<.00002`,
with zero cross covariance. First prove `P(t)<=U` throughout the first partial magnetic cell after handoff: bound
the prediction of `diag(.087^2,1e-6)` by its endpoint diagonals and absolute
cross term at `.04`, and check the positive residual diagonals and determinant.
The first magnetic correction can only decrease this covariance. Thus no
coincidence of the handoff and magnetic clocks is assumed. Thereafter,

`U - Riccati_.04(U) > 0`

by exact rational LDL.

Monotonicity of the regular scalar-observation Riccati map proves this bound
at every subsequent magnetic posterior. For every `0<=t<=.04`, bound the
two prediction diagonals and the absolute cross term by their endpoint
values. The determinant of the resulting lower comparison residual exceeds
`9.9983e-9`; both diagonals are positive. Therefore `P(t)<=U_prefix` at
**every** operation prefix, including roots between magnetic callbacks.

Every one-second window contains 25 consecutive magnetic corrections, spaced
`.04` seconds. Let O have normalized rows `(1,.02*t_i)`. A common time shift
changes neither its Gram determinant nor rank. Using `0<=t_i<=1.04` safely
covers arbitrary window phase and a preceding IMU root. Exactly,

`det(O^T O)=.02^2*25^2*(25^2-1)*.04^2/12 = .0208`,
`tr(O^T O) <= 25*(1+.02^2*1.04^2)`.

For the auxiliary Gaussian covariance model, the joint observation-noise
covariance V, excluding root uncertainty, includes **all** correlated process
noise. It obeys

`V <= v I`, `v=r+25*(sg^2*1.04+qb*1.04^3/3)`.

Indeed its process part is PSD and its trace is bounded by the sum of the
25 unconditioned heading variances. No process correlations are discarded.
Thus `J=O^T V^-1 O >= j I`, where `j=det/(tr*v)`.

Here the root-information calculation is legitimate: the triangular
innovations transformation whitens the same joint observation covariance
`O Pbar O^T+V` with the **actual** Kalman gains and innovation variances.
Its deterministic response to a root perturbation is exactly the complete
preceding corrected transport. Therefore the magnetic loss restricted to this
invariant pair is

`D_mag = O^T (O Pbar O^T+V)^-1 O = (J^-1+Pbar)^-1`.

This is information about the initial root; no endpoint batch-covariance
comparison is used. It is also checked independently against exact rational
sequential innovation losses with interleaved process noise.
Since `Pbar<=.075 I` in scales `(1,.02)`,

`D_mag >= (1/j+.075)^-1 I > 3 I`.

The exact lower bound is about 3.72037. It holds for every tail window,
with the carried covariance, actual gains and actual innovation covariance.
If the heading/axial-bias injection is instead taken about the physical body
down axis `(0,+/-4/5,3/5)`, the independent z group contributes `9/25` of
this matrix and the y group contributes a PSD matrix. Hence the service
floor is still greater than `27/25 > 1` (exact bound about 1.33933).
The required `T_M=1`, `mu_M=1` is therefore met under either down-axis
convention. No absence of magnetic service excuses the capture failure.

## 4. Capture contradiction and its precise scope

The literal common estimate is level. In either physical history,

`tilt_error(t)=acos(3/5) > pi/30`

for all times, including after the finite A21 release. No finite capture time
into the declared six-degree physical tilt region exists. More generally,
the two true down axes are separated by `2*acos(3/5)`. A deterministic
estimator receiving their identical sample histories cannot eventually stay
within any smaller common tilt radius than `acos(3/5)` for both histories.
This information obstruction is independent of arithmetic precision.
The source-uniform service certificate here is real arithmetic; no infinite
float32 covariance/clock-totality claim is substituted for it.

This **refutes the stated universal capture obligation** under the current
point-sample assumptions. It does not refute conditional local A21 stability,
prove covariance divergence, or show that a practical bound of arbitrary size
is impossible. The nominal filter can be contractive around its own level
trajectory while the admitted physical truth remains far from it. The prior
stationary detectability and nuisance covariance results remain valid.

## 5. Concrete missing premise; not silently adopted

Deleting exactly quiet motion does not remove these two nonzero histories.
A quantitative connection between sampled acceleration and physical velocity
increments is needed before continuous marine balance can prove sampled
attitude diversity or capture. For example, on every IMU cell require

`||h_k*a(t_k) - integral_[t_k,t_(k+1)] a(t) dt|| <= eta*h_k`.

Then, for any window made of complete cells and a constant world field,

`(1/T)||sum h_k*(a(t_k)-g) x B||`
`>= g*B_horizontal_min - (2*V_max/T+eta)*B_max`.

This follows directly by summing the cell errors and using
`integral a = Delta v`. A physical bound `||a'||<=J` suffices with
`eta=J*h_max/2`. As a concrete **unadopted** example, `J=100 m/s^3` gives
`eta=.3 m/s^2`; the displayed diversity lower bound at T=8 s is 21.47475,
strictly positive. The counterexample has cell error per unit time equal
to `2g/sqrt(5)`, and is excluded by this condition.

This repair still needs physical/sensor qualification and transfer from
physical force to the retained shipping trajectory; it does not itself close
AG loss, nonlinear retention, or the complete theorem. This premise would
narrow the current physical-history contract and needs physical justification.
Neither it nor its illustrative numerical value has been inserted into the
admitted domain or claimed as qualified.
