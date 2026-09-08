# Three complete-word storage routes

All three experiments use the same passively attached physical source capture:
H18 at source time 32.995–35.995 s and A21 at 1142.995–1145.995 s, each with
600 predictions and 600 accepted accelerometer corrections. They retain the
full 21-state finite map, actual covariance and anisotropic R_S, finite
quaternion resets, bias projection and common forcing history. Coefficients
are frozen at this one history. No replay optimization selects the metric.

`three-storage-routes.json` contains input and producer SHA-256 values,
maximizing physical root directions, 80-digit endpoint/prefix propagation,
operation energy accounting, information spectra and all failed metrics.
High precision rechecks the evaluated binary64 coefficients; it is not an
outward enclosure or proof for the nonlinear source family.

| Complete-word test | H18 | A21 |
|---|---:|---:|
| Centered motion endpoint rho | .999572659 | .932948382 |
| Centered motion maximum prefix rho | 1.003245112 | 1.000001212 |
| Additional full active-error endpoint rho | — | .997124968 |
| SI diagonal motion endpoint rho | 1072.993 | 6659.686 |
| Gravity/3 s diagonal motion endpoint rho | 11.183777 | 6.576082 |
| Full signed transported-information margin | .000427341 | .067051618 |
| Separate scalar information/remainder lower bound | −.380703 | −.257014 |

## Source-centered storage and separate budgets

The exact frozen factorization is `x_j = T_j x_0 + r_j`, where `r_0=0` and
`r_next=A*r+B*u`. Every physical forcing port uses the same template amplitude
alpha. The finite reset and projection coefficients remain in A and B.
The metric is the principal 18×18 block of the inverse **full** covariance.
For both modes the complete budget is

`sqrt(W_j) <= sqrt(rho_j)*sqrt(W0_motion) + sqrt(B_j)*|b0| + sqrt(F_j)*|alpha|`.

B measures the full map from all three initial bias coordinates, in m/s²;
F is the particular-response energy. No bias coordinate is deleted from the
recursion. The homogeneous motion ratio alone sets initial bias deviation to
zero; the separate B term restores arbitrary initial bias in this coefficient
inequality. The additional A21 full-error metric is reported separately and
cannot replace the bounded-bias motion theorem.

| Remaining point budget | H18 | A21 |
|---|---:|---:|
| Maximum squared initial-bias gain B | 3160.343 | 314671.332 |
| Maximum particular-response energy F | 4.13762 | 16.50404 |
| Sufficient .4 bias ball + template prefix bound | 600.633 | 51794.182 |
| Minimum point 30° chart storage level | 93.686 | 22824.342 |
| Sufficient bound / chart level | 6.41 | 2.27 |

The .4 bound here uses the capture's zero true bias and the unchanged closed
bias-estimate ball. It does not qualify other true-bias histories. These
sufficient bounds exceed the chart even before adding initial motion energy.
That makes this separated budget unusable for chart retention; it does not
show an admissible nonlinear trajectory leaves the chart.

The H18 limiting endpoint direction couples Y velocity, Y wave acceleration
and Y displacement. A21 is dominated by integrated Z displacement. Endpoint
and worst-prefix directions differ and are both retained. H18's limiting
prefix is the accelerometer projection at sample 6602; A21's is the wave
acceleration covariance floor at sample 228605.

## Physical scaling

The first diagonal storage uses one SI unit per coordinate. The second uses
fixed units `[1 rad, 1/3 rad/s, 3g m/s, 9g m, 27g m s, g m/s², g m/s²]`,
repeated on each axis. Three seconds is the declared complete-word duration;
g is the capture's gravity constant. Neither is fitted to observed extrema.
The full observer and covariance are unchanged when evaluating these storages.

Both physical metrics fail by large margins. A separate control conjugates
the transition and transforms the original information metric consistently;
its rho is invariant. Changing numerical coordinates cannot change a failed
mathematical inequality. These two diagonal metrics are rejected without
interval refinement. This does not reject every possible physical storage.

## Transported vector information

For every executed measurement, the exact finite homogeneous residual H is
transported by the full corrected transition **before** that measurement:
`J = sum T_before^T H^T R^-1 H T_before`. Actual accelerometer, magnetometer
and S-zero covariances are retained, with separate channel contributions.
Every prefix consumes only its own past measurements. The eta6 SI minimum
eigenvalues are 972.647 and 862.261 at this point.

The controlling identity is
`M0 - T_end^T M_end T_end = J + signed_remainder`.
It reproduces the endpoint margins in the table. The remainder explicitly
retains prediction, covariance floors, finite correction, reset and projection
effects; it is not assumed positive. Scalarizing J and the remainder
separately loses the favorable matrix alignment and produces negative lower
bounds. That scalar tactic fails; the full signed matrix remains feasible
at the captured point. This ledger restates the same complete-word margin,
so it is not an independent stability proof.

## Reproduction and status

Run the passive capture and attachment commands in
`.github/workflows/ou3-p4-physical-finite-map-feasibility.yml`, then:

```sh
OPENBLAS_NUM_THREADS=1 python3 tests/kalman_ou_iii/ou3_p4_storage_routes.py \
  --prefix /tmp/ou3-source-events \
  --attachment /tmp/ou3-connected-motion.json \
  --output /tmp/ou3-three-storage-routes.json
```

The workflow runs these three routes after the explicit BIAS2 endpoint/prefix
tests and uploads their output. Focused tests verify coordinate invariance,
the signed information identity, the separate bias budget and prefix causality.

Physical source admission, uniform BIAS2 separation, a uniform bounded
particular response and a nonlinear coefficient enclosure remain open.
Conditional canonical P3 remains at 1e-18 with H18/A21. All promotion flags
remain false; P5 cannot start. Current failure analysis and alternatives are
in `docs/ou3-proof-research-state.md`.
