# Acceleration regularity, sampled information and physical prediction supply

These results enter the existing finite-error proof in two places: physical
prediction mismatch contributes to the disturbance supply in
`V_next <= rho V + c_d ||d||^2`; sampled vector information is an input to
observability. Neither replaces the complete 21-coordinate corrected loss.

## Adopted physical condition

MARINE MOTION now requires acceleration to be locally absolutely continuous,
with `||a_dot(t)|| <= J_max = 100 m/s^3` almost everywhere. Equivalently it is
100-Lipschitz on every compact time interval. This concerns the same physical
acceleration in `p'=v, v'=a`; reference acceleration and mounting forces may
not be hidden by changing coordinates. Existing declared disturbances remain
charged separately. The displacement primitive continues without resets.
Quiet water remains admissible. There is no minimum wave amplitude.

The limit caps acceleration change at 0.6 m/s^2 across the longest admitted
0.006 s IMU interval. It is a selected certification limit, not a universal
law for all vessels or an assertion that every installation is qualified.
For physical scale, the primary reference generator
`oceanography-waves-lib` tag v1.2.1, commit
`e442150682f560384be427df4cc7815956a091c5`, represents CG translation as a finite
sum of harmonics with coefficients `c_n` and frequencies `omega_n`. Its
all-time triangle envelope is `sum omega_n^3 ||c_n||` (complex Euclidean norm).
Evaluating the unmodified `VesselRao::transfer`, `incidentHarmonics` and
`regularWaveHarmonics` for its 20 declared scenarios gives a largest envelope
about 29.1126 m/s^3, for the PM 11.4-s case. This floating coefficient screening
supports the scale of the selected limit; the proof below does not depend on
that screening or claim measured-hull qualification. Actual deployments must
certify their continuous-history jerk bound, not infer it from finite samples.

The old 200-Hz witness has jerk amplitude `400*pi*2g/sqrt(5) > 10500 m/s^3`.
It is excluded. Its exact obstruction remains useful for showing why removing
this condition would invalidate the six-degree capture claim.

## 1. Sharp sampling bound on nonuniform cells

For a cell of length h, integration by parts gives

`h(a(0)+a(h))/2 - integral_0^h a(s) ds = integral_0^h (s-h/2) a_dot(s) ds`.

Therefore the norm is at most `J h^2/4`. This is sharp: a triangular tent
with slopes +J and -J attains equality. In contrast, left-endpoint sampling
has error at most `J h^2/2`. Absolute continuity is essential to this proof;
a bound on an a.e. derivative alone would permit jumps or singular variation.

For complete nonuniform cells of total duration T, define endpoint weights
`alpha_0=h_0/(2T)`, `alpha_N=h_(N-1)/(2T)`,
`alpha_i=(h_(i-1)+h_i)/(2T)`. They are positive and sum to one. Since
`integral a = v_end-v_start`,

`||sum alpha_i a(t_i)|| <= 2 V_max/T + J sum h_i^2/(4T)`

`<= 2 V_max/T + J h_max/4 = 11/T + .15`.

This weights existing samples for the proof; it does not substitute a
trapezoidal integrator for any shipping operation. It also proves the
necessary quadrature tolerances used by the finite trace audit. A finite
trace that passes remains insufficient for all-time membership.

## 2. The entire fixed-attitude stationary-sample ambiguity is excluded

Let true body-to-world R be constant and suppose every accelerometer sample
on a complete-cell window equals `-g e_z`, including physical bias and sensor
residual. The physical sensor equation implies

`a_i = g(e_z-R e_z) - R(b_a,i+n_a,i)`.

Taking the weighted mean and using the preceding bound gives, with true tilt
`theta=angle(e_z,R e_z)`,

`2g sin(theta/2) <= 11/T + .15 + B_a + N_a`.

At T >= 32 s the right side is at most 1.018916605 m/s^2. The exact rational
certificate bounds `2g sin(3 degrees)` strictly above this using
`pi>3.14159` and `sin(x)>=x-x^3/6`. Consequently `theta<6 degrees`.
This excludes all constant-attitude stationary-looking aliases of at least
six degrees over such windows, with the full allowed accelerometer residuals,
not just the original sinusoid. It does not prove capture for moving R or for
nonstationary records of the actual adaptive observer.

## 3. Joint physical vector information, including finite rotation errors

For a constant world magnetic field B with `||B||<=75` and
`||B cross e_z||>=15`, transport the measured accelerometer samples into
world coordinates using true R and divide by g; call them u_i. Transport
actually applied magnetic samples by true R and divide by 75; call them b_j.
The use of true R defines an analytical coordinate, not a quantity supplied
to the estimator. Give accelerometer rows the weights above, and applied
magnetic rows any nonnegative weights beta summing to one. Magnetic service
ensures these events exist; no information is inferred from callback cadence.
Charge both the 5-uT hard-iron and 2-uT measurement residual envelopes.

Writing weighted means u_bar,b_bar and e=e_z,

`||u_bar+e|| <= epsilon = (11/T+.15+B_a+N_a)/g`,

`||b_bar-B/75|| <= r = 7/75`.

The full three-dimensional matrix (each skew factor is rank two)

`G = sum alpha_i [u_i]_x' [u_i]_x + sum beta_j [b_j]_x' [b_j]_x`

obeys `G >= [u_bar]_x'[u_bar]_x+[b_bar]_x'[b_bar]_x` by variance positivity.
No attitude coordinate or cross term is removed. For any two vectors u,b,
the two smaller eigenvalues of their joint matrix have sum
`||u||^2+||b||^2` and product `||u cross b||^2`; the remaining eigenvalue is
that sum. Thus its least eigenvalue is at least product/sum.

The triangle inequality gives

`||u_bar cross b_bar|| >= 1/5 - epsilon - r(1+epsilon)`.

For T >= 64 s this is positive, and consequently

`G >= gamma I_3`,

`gamma = (1/5-epsilon-r(1+epsilon))^2 / ((1+epsilon)^2+(1+r)^2)`

`>= 0.0000629714509407403` (the artifact retains the exact rational value).

This also supplies a finite-angle identity. For a common rotation Q of angle
phi and unit axis z, Rodrigues' formula gives

`sum alpha_i ||(Q-I)u_i||^2 + sum beta_j ||(Q-I)b_j||^2`

`= 4 sin(phi/2)^2 z' G z >= 4 gamma sin(phi/2)^2`.

It does not require a first-order attitude approximation. For different
rotations Q_i at different rows, stack the weighted residuals. The reverse
triangle inequality bounds their norm below by
`2 sqrt(gamma) |sin(phi_root/2)|` minus the norm of the stacked variations
`(Q_i-Q_root)u_i` and `(Q_j-Q_root)b_j`. Those variations and nuisance-state
cancellation must still be bounded using the actual observer dynamics.
In particular this physical Gram matrix is not the nominal, innovation-
weighted, corrected full-state `D_word`; equating them would repeat the
previous information-lifting error.

## 4. Physical mismatch in the literal OU prediction storage

Let z=(v,p,S,a) denote physical LIN coordinates with the same once-anchored S.
On cell k freeze the actual realized `lambda=1/tau_k`. Without any OU law on
truth, the physical identity is

`a_dot = -lambda a + u`, with `u=a_dot+lambda a`.

The ideal integrated-OU transition F_ideal therefore gives exactly
`z_next=F_ideal z + integral exp(A(h-s)) E_a u(s) ds`.
Its process covariance uses `W=2 lambda Sigma_aw`. The endpoint minimum-action
inequality, proved by weighted Cauchy-Schwarz, yields

`d_ideal' Q_ideal^-1 d_ideal <= integral u' W^-1 u ds`.

All four LIN blocks, three spatial axes and process cross covariance remain
in this inequality. For the source small-x branch, put
`d_source=d_ideal+(F_ideal-F_source)z`. Only physical a appears in the latter
term. The existing source kernel certificate gives
`Q_source >= (1-epsilon_Q) Q_ideal` and a coefficient-action bound
`||(F_ideal-F_source)z||_(Q_ideal^-1)^2 <= A_max^2 c_F`.
The closed-form branch has zero transition defect. Young's inequality with
`eta=1/100` gives the uniform shipping bound

`d_source' Q_source^-1 d_source`

`<= ((1+eta) integral u'W^-1 u ds + (1+1/eta) A_max^2 c_F)/(1-epsilon_Q)`.

The committed certificate also supplies a finite, deliberately conservative
box bound using `Sigma_aw >= .05^2 I`, `lambda in [1/12,50]` and
`||u|| <= J+lambda A`. The maximum of `(J+lambda A)^2/lambda` occurs at an
endpoint because it is convex. Summing the resulting per-cell action uses
`h>=.004`; no independent resampling of a is allowed.

This enters the **same full covariance metric**. For the full prediction
`P_minus=F P F'+Q_total`, with the literal LIN process block in Q_total and
any actual PSD sync retained, the variational identity gives

`(F e+E_LIN d_source)' P_minus^-1 (F e+E_LIN d_source)`

`<= e'P^-1 e + d_source'Q_source^-1 d_source`.

One proof regards `P_minus` as the covariance of the sum of independent
auxiliary root and process vectors and minimizes their total quadratic
energy subject to that sum. This is an algebraic covariance comparison, not
a stochastic assumption on physical truth. It preserves root cross covariance
and does not replace the source Joseph arithmetic. Other physical forcing,
nonlinear corrections, projection, reset and float32 supply still have to be
composed with a verified strict complete-word loss and prefix retention.
