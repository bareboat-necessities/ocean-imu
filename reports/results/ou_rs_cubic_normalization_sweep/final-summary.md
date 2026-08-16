# Cubic effective-rS normalization sweep

This experiment fixes the candidate adaptation law to

\[
r_{S,\mathrm{eff}} = K\,\sigma_{aw}\,\tau^3
\]

and varies only the normalization `K`. The full OU-III Adaptive suite is used: four stationary JONSWAP seas, four stationary PM-Stokes seas, the standard non-stationary transition, and the same ten predeclared paired seed triplets. Each arm contains 90 paired validation rows.

The nominal reference point is `tau=2.17904091 s`, `sigma_aw=0.724445343 m/s^2`, so `rS_nom = 7.49552016 K`.

## Results versus deployed law

Negative percentages are improvements.

| nominal effective rS [m s] | K | stationary 3-D | stationary Z | blend 3-D | blend Z | whole-transition 3-D | whole-transition Z |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.1606 | 0.154836 | -4.56% | +1.54% | +3.62% | +8.52% | -0.22% | +4.26% |
| 1.3000 | 0.173437 | -2.56% | +0.51% | +1.92% | +4.50% | +0.16% | +1.82% |
| 1.4000 | 0.186778 | -0.80% | +0.22% | +1.26% | +2.25% | +0.90% | +0.65% |
| 1.5000 | 0.200120 | +1.18% | +0.22% | +0.97% | +0.42% | +1.93% | -0.15% |
| 1.6000 | 0.213461 | +3.33% | +0.46% | +1.00% | -1.05% | +3.22% | -0.62% |
| 1.7000 | 0.226802 | +5.61% | +0.90% | +1.29% | -2.22% | +4.71% | -0.82% |
| 1.8600 | 0.248148 | +9.48% | +1.95% | +2.20% | -3.55% | +7.42% | -0.68% |

The repeated `rS_nom=1.160579` arm reproduces the prior cubic experiment exactly: paired values for X/Y/Z/3-D displacement, transition segments, tau and sigma are bit-for-bit identical.

## Decisive finding

A local quadratic fit through the middle normalization points gives:

- stationary 3-D crosses deployed performance at approximately `rS_nom = 1.441 m s` (`K = 0.19223`); tighter values improve stationary 3-D, looser values worsen it;
- transition-blend Z crosses deployed performance at approximately `rS_nom = 1.526 m s` (`K = 0.20355`); values below this worsen blend Z;
- transition-blend 3-D has its fitted minimum near `rS_nom = 1.548 m s` (`K = 0.20657`), but even there remains approximately `+0.93%` worse than deployed.

Therefore there is no constant cubic normalization that both preserves the stationary 3-D advantage and avoids the transition penalty. The required stationary-3-D and transition-Z regions do not overlap, and transition-blend 3-D remains worse throughout the tested cubic family.

## Practical interpretation

- If only stationary full-3D performance matters, the tighter cubic region is beneficial. `rS_nom=1.30` is a useful compromise: stationary 3-D improves 2.56% while stationary Z RMS worsens only 0.51%; all eight stationary sea-state mean 3-D results improve, with 79/80 individual sea/seed pairs improving.
- If a single constant cubic law is forced across stationary and transition operation, about `rS_nom=1.50` is the least-bad balance among the tested points, but it still worsens stationary 3-D by 1.18% and transition-blend 3-D by 0.97%, so it does not justify replacing the deployed law.
- The data support a two-regime design rather than another global normalization: use the tighter near-cubic law in settled/stationary conditions and relax/revert the regularizer during detected sea-state transitions. That transient scheduler requires its own paired validation before production use.
