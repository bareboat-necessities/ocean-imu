# Finite-error carried-storage feasibility and the nonlinear bridge

## Question and result

Does the prescribed storage survive actual nonzero errors, reset/projection,
and held-to-active bias release? The native finite-history experiment says yes
on the tested stationary BIAS0 histories. It does not establish the supremum
on the nonlinear BRMM domain.

`finite_error_storage_diagnostic.py` runs the unmodified shipping wrapper from
construction. The physical source has identity attitude, zero translation and
physical biases, gravity-only accelerometer input and a horizontal 32 uT field.
After actual Live, magnetic-reference acquisition and internal bias unlock,
a 100-sample gyro disturbance creates the error through normal inputs. The
predeclared amplitudes are 0.001 and 0.01 rad/s, each along each body axis;
both lie below the declared BMI270/ICM gyro residual caps. The disturbance
then ends. Magnetometer callbacks continue every eight global 5 ms samples.
This describes attempted calls; accepted informative recurrence is not inferred.

Each of H18, A21 and actual H18-to-A21 release runs for 1,800 additional samples
without reseeding; release is before sample 601. There are 18 histories and
54 three-second words. No estimator state, covariance, gain or metric is fitted
or installed. Shipping prediction, Joseph arithmetic, reset, projection,
frontend and tuner execute normally. The probe checks Live, dormant vibration
guard, finite state, disabled lever arm, unchanged magnetic reference and
absence of magnetic refinement at each measured sample boundary.

The actual error uses the world-to-body convention:

`c = 2 vec(q_true_WB conjugate(q_hat_WB)) / scalar`,

followed by true-minus-estimated gyro bias, v, p, centred S, a_w and accelerometer
bias; the last three joint coordinates are true physical beta, here zero.
The complete carried P21, including held-bias cross terms, defines the energy.

The worst measured finite ratio is **0.9744158786871832**, for the z-axis
0.001 rad/s prefix, shared by all three modes. This is a sampled-history ratio,
not a maximizing direction or a bound comparable as a supremum to 0.97715855.
All 72 word endpoints are exactly finite binary32 SPD by rational LDL.
All 32,418 measured sample boundaries pass floating Cholesky, with no symmetry
repair or eigenvalue floor. Their largest prefix/start storage ratio is
1.0000000000000013 (numerical evaluation precision). Internal intermediate
writes and startup/pulse covariance boundaries are not covered by that audit.

The report `reports/results/ou3_alt_storage/finite-errors.json` retains every
endpoint error, exact minimum pivot, 60-digit endpoint energy/ratio, trace
hash, compiler and shipping source-tree identity. The 60-digit calculation
re-evaluates serialized binary64 Cayley coordinates: it is not an enclosure
of their division or of target-machine execution. Nonzero physical motion,
BIAS1/2, arbitrary continuations and informative-acceptance bounds remain open.
Zero sensor disturbance after the prefix does not eliminate machine supply.

## Conditional nonlinear comparison lemma

This gives a quantitative target for the remaining proof, not measured constants.
Fix a history and a declared word domain. Let N0,N1 be positive reference
metrics and A a linear comparison operator. The input e denotes the declared
unsupplied coordinates; the output retains the full storage coordinates. True
physical beta and any held-bias component assigned to supply must enter D and
the remainder explicitly. A may therefore be rectangular. More generally the
input premise below is V_now >= (1-di)||e||_N0^2, which must include any
initial cross terms rather than discarding them. Suppose **uniformly on that domain**:

- `||A e||_N1 <= a ||e||_N0`;
- the full finite word satisfies `F(e,d)=A e+r`, with
  `||r||_N1 <= eps ||e||_N0 + b D`;
- the actual input energy satisfies `V_now>=(1-di)||e||_N0^2` and
  the output metric satisfies `M1(e,d)<=(1+do)N1`, where `0<=di<1`,
  `do>=0` (for a full input vector, `M0>=(1-di)N0` suffices);
- D includes all physical/model/rounding supply used in the remainder bound.

These hypotheses concern the same carried source/covariance history, including
all branch changes and release. A numerical tangent on one quiet history does
not establish them. Choose eta>0. Triangle inequality and weighted Young give

`V_next <= rho V_now + C D^2`,

`rho = (1+do)/(1-di) (1+eta) (a+eps)^2`,

`C = (1+do)(1+1/eta) b^2`.

Proof: bound the output norm by `(a+eps)||e||_N0+bD`, square with
`(x+y)^2<=(1+eta)x^2+(1+1/eta)y^2`, multiply by `(1+do)`, and use
`||e||_N0^2<=V_now/(1-di)`. No differentiability of the rounded map is assumed;
its quantization defects must be included in the remainder/supply bound.
This also allows a finite offset at zero error; it cannot be silently set to zero.

An exact rational *illustration*, a=989/1000, eps=1/1000,
di=do=1/500 and eta=1/200, gives rho<1. The a value is above the square root
of the quiet tangent measurement on its unsupplied motion subspace; it is not
a certified operator bound on either the full joint24 state or a nonlinear
source domain. `nonlinear_budget()` computes these coefficients exactly and rejects
invalid metric budgets. The example explains the available scale: roughly
0.2% metric variations and a 0.1% norm remainder would fit inside the observed
margin. None of those hypotheses has been proved for shipping histories.

## Controlling next experiment

Generate nonzero BRMM motion and physical BIAS1/2 histories through the same
construction, compute the actual physical truth and supply, and evaluate the
whole-word residual in the prescribed metric. Retain magnetic acceptance and
all covariance/branch changes. The tested pulse directions do not cover the
limiting translation direction of the earlier tangent supremum. A finite ratio
above one would require checking the supply term before declaring a theorem
failure. Source-uniform remainder/metric comparison and indefinite finite-SPD
preservation remain the central blockers; no P4/P5 gate is promoted.

Reproduce from the repository root:

```sh
EIGEN_INCLUDE_DIR=/path/to/eigen OPENBLAS_NUM_THREADS=1 PYTHONPATH=. \
python3 tools/stability/ou3_alt_contraction/finite_error_storage_diagnostic.py \
  --work-directory /tmp/ou3-finite-errors \
  --output reports/results/ou3_alt_storage/finite-errors.json
```
