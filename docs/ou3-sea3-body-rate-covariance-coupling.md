# SEA3 body-rate covariance coupling

The #482 translational RAO theorem intentionally left the parent vessel's rotational response under the separate Normal-Live body-rate source bound. `tools/ou3_sea3_body_rate_covariance_coupling.py` supplies the missing analytical coupling interface without declaring or fitting one universal rotational hull model.

For a candidate rotational displacement RAO

```text
||r(f,theta)||_2 <= K min(1,(f_c/f)^q),  q >= 1,
```

where `r` maps sea-surface elevation to body attitude in rad/m, body rate satisfies

```text
omega_hat = i 2 pi f r eta_hat.
```

For every frequency and heading,

```text
(2 pi f)^2 ||r||^2 <= K^2 (2 pi f_c)^2.
```

Since `m0 = Hs^2/16`, the rotational response has

```text
tr Cov[omega_rad/s] <= Hs^2 K^2 (2 pi f_c)^2 / 16.
```

The radian-to-degree conversion cancels pi exactly, leaving the particularly simple theorem quantity

```text
tr Cov[omega_deg/s] <= 8100 Hs^2 K^2 f_c^2.
```

No spectral quadrature or JONSWAP-shape constant is needed. Conditional on the rotational envelope, the bound covers every admitted JONSWAP gamma, arbitrary directional spreading, arbitrary complex response phase, and PSD-consistent cross-axis coupling.

For the 20-minute / 200 Hz horizon (`N=240000`) the existing finite-horizon concentration certificate allocates half of the 5% stochastic budget to the 30 deg/s body-rate event. The resulting required body-rate covariance trace is just below `8.3333333333 deg^2/s^2` because the proof threshold is deliberately rounded downward.

The new certificate keeps the coupling test in squared form,

```text
(Hs K f_c)^2 <= v_rate / 8100,
```

so the PASS decision needs no square root. As a nonphysical interface witness only, `Hs=8.5 m`, `K=1/8 rad/m`, `f_c=0.03 Hz`, `q=1` gives an outward body-rate covariance bound near `8.23 deg^2/s^2`, inside the 20-minute threshold. The same `K` at `f_c=1.2 Hz` fails by a large margin. This demonstrates that rotational coupling, like translation, must keep sea severity and response bandwidth correlated rather than independently maximizing them.

The `K=1/8` witness is **not** a claimed vessel RAO and does not qualify a physical hull population. The remaining SEA0 obligation is to attach a validated rotational response envelope (or a tighter direct body-rate covariance enclosure) for the intended vessel class, then compose that result with the already-certified acceleration covariance candidate in `ou3_sea3_finite_horizon_concentration.py`.

No filter constant, P1 cap, P2 source interface, P3 gate, P4 gate, or theorem promotion flag is changed by this certificate.
