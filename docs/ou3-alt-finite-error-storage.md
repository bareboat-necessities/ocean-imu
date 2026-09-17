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

## Coupled moving-source result

`driven_storage_diagnostic.py` now tests a nonzero physical source from actual
construction, independently for H18, A21 and actual release, BIAS0/1/2, and
omega=1/2 or 1 rad/s. The single fixed source is

`p(t)=(1/4)cos(omega t) d`, `d=(0,3/5,4/5)`,

with its exact derivatives v,a, identity physical attitude, zero physical gyro
bias, and constant horizontal magnetic field. The physical primitive is
`U(t)=sin(omega t)d/(4 omega)`. The actual Live transition sets the one origin:
`S(t)=U(t)-U(t_Live)`; the measurement root and later words never reset it.

Each bias history acts along body y from t=0:

| Family | Physical beta_y(t), SI units | Persistent family parameter |
| --- | --- | --- |
| BIAS0 | 0.02 + 0.005 sin(2 pi t/600) | zero GM channel, tau=600; turn-on offset plus non-GM channel |
| BIAS1 | 0.08 exp(-t/1200) + 0.015 sin(2 pi t/600) | one root, tau=1200, period=600 |
| BIAS2 | 0.05 + 0.01 sin(2 pi t/600) | non-relaxing phi_true=1 |

The analytic member certificates consume the existing physical primitive and
separate bias contracts. With 3<pi<22/7 they prove all-time magnitude, derivative
and increment bounds, rather than checking only replay samples. For BIAS0,
`|w|<=h(rate+B/600)`; for BIAS1 the exponential root cancels and
`|w|<=h*0.015*(2pi/600+1/1200)`; for BIAS2 `|w|<=h*rate` with phi=1.
The native h is the binary32 value of 0.005 and is less than 1/200, used in
these bounds. The wave certificate gives `|S|<=1/(2 omega)` for any Live time.
The frequency range is inside [0.018,0.88] Hz. This is a sufficient member
construction, never a restriction of COMPLETE-BRMM to harmonics.

All 18 histories reach the unlocked diagnostic root (samples 8,360..9,344).
They produce another 54 complete words and 72 exactly SPD endpoints. The
maximum **raw** ratio is **4.104252247690883**, for A21/BIAS0/omega=0.5 in
samples 600..1200. Therefore the quiet unforced ratio does not extend to these
forced words. This is not a counterexample to dissipativity with supply.

To expose a necessary gain without fitting a certificate, predeclare a candidate
physical supply vector `xi=(p,v,S,a,beta)`, each coordinate divided by one in
its stated SI unit. Its analytic, all-time squared envelope is

`D_phys^2 = (1+omega^2+omega^4)/16 + 1/(4 omega^2) + B^2`,

where B=1/40,19/200,3/50 for the respective bias members. Proof: apply the
unit-direction amplitude bounds to p,v,a, the endpoint-difference bound to S,
and the all-time magnitude bound to beta; add their squared bounds. This
supply contains physical quantities only, not error, P, or a replay-fitted
state radius. The choice of unit scales is explicit; C depends on it.

For this candidate convention, every inequality `V1<=rho V0+C D_phys^2` must
satisfy `C>=max(0,V1-rho V0)/D_phys^2`. The report evaluates that **necessary**
condition at predeclared rho targets 0.95,0.98,0.99; it does not select C to
claim success. At rho=0.98, the strongest observed requirements are:

| Family | Largest raw ratio | Necessary C, approximately |
| --- | ---: | ---: |
| BIAS0 | 4.104252 | 390.7913 |
| BIAS1 | 1.658480 | 417.3145 |
| BIAS2 | 2.209533 | 404.3313 |

The largest required C occurs in H18/BIAS1/omega=1, samples 1200..1800:
V rises from 511.66465 to 687.77272. The signed allocation `e_i(P^-1 e)_i`
places about 442..444 in the held accelerometer-bias group, while the position
allocation rises from 62.87 to 216.73 and S from -8.69 to 37.71. These are
signed allocations retaining cross terms, not independent positive energies.
The report includes every endpoint allocation. Dropping held-bias output
energy would materially change the tested storage and is not permitted.

These are host diagnostic values, not outward-rounded certificates. They
include all 24 output coordinates, including beta energy and the complete P21
inverse. The largest raw ratio and the largest C need not be the same word.
If a theorem introduces a separate additive machine budget, that budget must
be subtracted before claiming this same lower bound on its physical C.

The analytic wave and bias bounds are proved for the ideal functions. Native
sin/exp, physical-state serialization and binary32 sensor conversion still
need rigorous attachment to the target arithmetic theorem. The trace reports
sensor conversion residuals; that observation is not a libm qualification.
Attempted magnetic callbacks do not certify informative accepted recurrence.
No universal rho, upper bound on C, nonlinear basin, indefinite SPD preservation
or P4/P5 gate is established by this diagnostic.

## Controlling next proof obligation

Bound the complete correlated nonlinear source-to-storage gain using the
canonical word's forcing channels and compare it with the necessary values
above under the same supply normalization. In particular attach the physical
OU mismatch, S=0 reference supply, bias mismatch, magnetic ancestry and machine
roundoff together. The original remainder/metric lemma remains valid, but its
b coefficient cannot be guessed from the quiet tangent margin. An alternative
source-driven reference theorem would additionally need a separate bound from
that reference to physical truth; it must not erase persistent tracking error.

Reproduce from the repository root:

```sh
EIGEN_INCLUDE_DIR=/path/to/eigen OPENBLAS_NUM_THREADS=1 PYTHONPATH=.:tools/stability \
python3 tools/stability/ou3_alt_contraction/driven_storage_diagnostic.py \
  --work-directory /tmp/ou3-driven-storage \
  --output reports/results/ou3_alt_storage/driven-words.json
```
