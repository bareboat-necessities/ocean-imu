# Conditional informative magnetic-service theorem

## Scope and selected hypothesis

The unconditional COMPLETE-BRMM source class is unchanged. The strict
full-heading conclusion applies only to a separately named Normal-Live regime
with informative magnetic service, after actual Live entry and heading
acquisition. This adds an explicit hypothesis to the ALT conditional theorem;
it does not assert that every admitted source supplies heading information.
The original P2/P3/P4/P5 route and P3=1e-18 are unchanged.

Declare positive constants T_B, L_B, h_B and alpha_B, with L_B >= 2 T_B.
For every overlapping interval [s,s+L_B] contained in this regime:

1. Use only actual accepted, gauged magnetic updates. Rejected callbacks,
   pending first-north acquisition and attempted service are insufficient.
2. Transport each actual whitened sensitivity through the full preceding
   same-history event product to the interval's initial heading/axial gyro-bias
   frame: B_i = R_i^(-1/2) H_i Phi_(s,i) E_s. The frame is orthonormal in
   the declared coordinates (radians and radians/second); arbitrary rescaling
   cannot manufacture an information floor. At finite error the corresponding
   incremental/nonlinear relation must be supplied, not a derivative of
   discontinuous floating-point rounding.
3. Informative events have ||B_i[:,0]|| >= h_B. Gaps between such events,
   including the initial and terminal boundary gaps, are at most T_B.
4. Their accumulated information satisfies sum B_i^T B_i >= alpha_B I_2.
   Counts or periodicity alone do not imply this rank-two inequality.

`InformativeServiceAssumption` makes these parameters and quantifiers explicit.
Its finite-window audit is a diagnostic, not proof of infinite recurrence or
stability. The service property may be assumed for a deployment or established
from additional source/runtime guarantees. It is no longer necessary to prove
it for every unconditional COMPLETE-BRMM history before stating this conditional
target. The whole-theorem gate still cannot pass from the assumption alone.

A nonzero horizontal magnetic field, regular sampling, bounded disturbances
and successful gates explain why this regime is physically relevant. They do
not guarantee it during saturation, interference, rejection or loss of heading.
The native experiment retains the declared horizontal field B=(32,0,0) uT and
existing physical input contracts. Its 40 ms, 200 ms and 1 s callback spacings are a
predeclared experiment grid, not new hardware guarantees. Only 40 ms retains
MAG-CALL-SCHEDULE-v1. The slower schedules and no-service outage are explicitly
classified as stress experiments outside that callback profile. They do not
qualify the declared theorem, even if the finite accepted-service audit passes.
The theorem retains the original attempted-call premise and adds informative
accepted service. Real rejection patterns at 40 ms attempted cadence require
nonzero-residual tangents/finite maps beyond the current probe. The numerical audit
instance T_B=1 s, h_B=0.01, alpha_B=1 and L_B=3 or 6 s is explicitly assumed;
these constants were chosen before measurement, not fitted to its extrema.
An infinite admissible history satisfying this instance is still an obligation.

## Optional history-wise target

`ou3-alt-nonuniform-capture-theorem.md` proves the conditional comparison result
when finite capture is assumed separately for each history and the strict tail
factor/coercivity constants may depend on that history. It needs no numerical
maximum capture time. This is a weaker eventual practical boundedness target;
it does not close shipping finite-error, coercivity or retention premises.
Merely having each time-varying factor below one is insufficient.

## Corrected stability target

Keep the prescribed compatible law M(P)=diag(P21^-1,I3). On every admitted
same-history service superword, the missing finite-increment certificate is

`V_(j+1) <= rho_B V_j + c_B ||d||²_[j,j+1],  rho_B < 1`,

with uniformly coercive endpoint storage and a bound on every intervening
prefix. The six canonical supplied coordinates, all output coordinates,
H18/A21/control edges, nonlinear correction, covariance/reset defects and
BIAS0/1/2 remain represented. Accepted service is one premise of this inequality;
it is not a proof of the inequality.

If the inequality and the prefix bound V(t)<=K V_j+c_p D² hold with uniformly
bounded supplies D, induction gives

`V_j <= rho_B^j V_0 + c_B D² (1-rho_B^j)/(1-rho_B)`.

Uniform coercivity then converts this to a state bound, and the prefix bound
extends it between superwords. Neither an endpoint spectrum nor growth of P
can substitute for these premises. No numerical rho_B is certified here.

## Ungauged regime

The heading gauge alone can be quotiented only if the actual map descends to
it. Test Q_after^T A N_before=0 before deleting coordinates. Even when this
passes on a quiet word, axial gyro bias remains neutral. It is not another
heading gauge and cannot simply be discarded or called a bounded supply.
The pair (heading, axial bias) has unipotent centre dynamics; centre-growth and
its coupling into the transverse system must remain in the theorem.

For comparison, the diagnostic also removes that pair algebraically and tests
invariance. This second calculation is a transverse/centre decomposition on
the measured word, not a physical quotient theorem. General moving histories
may fail its invariance condition. Both calculations report point leakage,
motion spectra and identity-storage ratios with all remaining output rows.
Identity-storage failure does not disprove the existence of compatible storage.
The compatible quotient metric is the minimum energy over the removed fibre:
`M_q = Q^T M Q - Q^T M N (N^T M N)^(-1) N^T M Q`.
Taking a principal submatrix would not represent that quotient norm.

The no-service native control begins at the same reached gauged root as the
serviced histories and then omits callbacks. It measures a finite magnetic
outage, not a new proof of pre-gauge startup or a universal no-magnetometer
history. The existing exact ungauged unipotent proof remains authoritative.

## Measured comparison

`reports/results/ou3_alt_storage/service-regimes.json` records twelve native
histories and 36 window measurements. Each entry below is the worst across
H18, A21 and the actual H18-to-A21 release, including both 3-second windows.
These are local tangent storage ratios, not certified universal bounds.

| Callback spacing | Runtime scope | 3-second rho | 6-second rho |
| --- | --- | ---: | ---: |
| 40 ms | Retained call profile | 0.9765963647 | 0.9509988673 |
| 200 ms | Outside-profile stress | 0.9768436096 | 0.9509988673 |
| 1 s | Outside-profile stress | 0.9798998419 | 0.9540771044 |
| No calls | Outside-profile outage control | 0.9996640645 | 0.9993281209 |

All nonzero schedules pass the finite informative-service audit; only 40 ms
also satisfies the retained callback profile. Actual pair information has a
minimum eigenvalue above 26,324 on the retained 3-second windows. Binary64
ratios and the independent 60-digit terminal checks differ by less than 3e-15.
All observed/plain sample states agree bit for bit.

The no-service control has zero point fibre leakage for both tested removals.
Removing only heading leaves motion spectral radius exactly 1. Removing the
heading/axial-bias centre pair gives worst transverse spectral radius
0.9739828093 and compatible transverse storage ratio 0.9766752765 over three
seconds; the six-second compatible ratio is 0.9509988673. This supports a
transverse-plus-centre architecture on this quiet history. It does not turn
axial bias into a physical gauge or prove invariance for moving sources.

The apparent full-storage contraction without any magnetic information is an
important negative control: P can grow along the unobserved pair and reduce
its inverse-covariance weight. The exact neutral/unipotent obstruction remains.
Uniform coercivity is indispensable, even when every measured ratio is below 1.

Next falsifiable work is the nonzero-motion, finite-error map with actual
rejected events at retained 40 ms cadence, and the quotient's transported
physical gauge/centre coupling on that same history. Uniform coercivity,
all BIAS families, every hybrid edge and source/target arithmetic remain open.
Do not start interval refinement of a stationary point and label it universal.

## Reproduction and evidence boundary

Run from the repository root with Eigen installed:

```sh
PYTHONPATH=.:tools/stability OPENBLAS_NUM_THREADS=1 python3 -m \
  tools.stability.ou3_alt_contraction.service_regime_diagnostic \
  --work /tmp/ou3-service-regimes \
  --output /tmp/ou3-service-regimes.json
```

The experiment compiles observed and plain native probes, drives both from
construction, and requires bit-identical sample states for every schedule.
It measures two consecutive 600-step windows and their composed 1200-step
superword from the same retained native trajectory. No intermediate covariance,
frontend, tuner or state is reseeded. Actual accepted H/R rows and preceding
updates supply the information audit. The report records maximizing vectors,
operation-wise energy changes and an independent 60-digit terminal check.

All measurements remain zero-residual local tangents on the stationary source.
They do not cover finite errors, general motion, all bias histories, arbitrary
word phases, infinite duration, or unrepresented magnetic reference refinement.
