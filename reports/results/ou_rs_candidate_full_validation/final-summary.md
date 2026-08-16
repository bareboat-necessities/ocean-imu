# Full paired validation of candidate effective-rS laws

This report compares the current deployed OU-III adaptive regularizer against two candidate laws on the complete ten-seed validation suite. Production filter source was not changed; the candidate law is injected only by the experiment harness.

## Protocol

- 4 stationary JONSWAP seas + 4 stationary PM-Stokes seas + the standard non-stationary JONSWAP transition.
- 10 identical predeclared wave-phase / IMU-noise / initialization seed triplets per scenario.
- 1200 s replays; trailing 900 s scored; transition interval 420--780 s.
- 90 paired rows per arm.
- Candidate effective laws:
  - cubic: rS_eff = 0.1548363522 sigma_aw tau^3
  - fitted: rS_eff = 0.1667018769 sigma_aw tau^2.9052
- Both candidates pass through the independently measured nominal full-3D optimum rS_eff = 1.160579 m s at tau = 2.17904091 s and sigma_aw = 0.724445343 m/s^2.
- The existing adaptive tau(t), sigma_aw(t), and pseudo-measurement cadence are unchanged.

The deployed, cubic, and fitted artifacts each contain exactly 90 unique paired rows. Pair keys match exactly. The maximum paired differences in tau_applied_s and sigma_applied_mps2 between deployed and either candidate are exactly zero. Simulator return-code / quality-gate patterns are also identical arm-to-arm, so the regularizer-law comparison does not alter which existing validation cases trip those gates.

## Eight stationary seas

Statistics below first average the eight stationary scenarios within each seed triplet and then compare the ten paired seed aggregates. Intervals are descriptive paired bootstrap 95% intervals over those ten seed aggregates.

| metric | deployed | cubic p=3 | cubic change | fitted p=2.9052 | fitted change |
|---|---:|---:|---:|---:|---:|
| Z error [%Hs] | 4.4887 | 4.5307 | +0.0420 pp [ +0.0108, +0.0696 ] | 4.5410 | +0.0522 pp [ +0.0193, +0.0816 ] |
| 3-D RMS [m] | 0.39519 | 0.37682 | -4.56% [ -5.40%, -3.83% ] | 0.37386 | -5.28% [ -6.34%, -4.35% ] |
| X RMS | -- | -- | -2.55% | -- | -2.81% |
| Y RMS | -- | -- | -8.58% | -- | -10.14% |
| Z RMS | -- | -- | +1.54% | -- | +2.06% |
| yaw RMS | -- | -- | +0.04% | -- | +0.04% |

The 3-D RMS improvement has the same sign on all 10 seed aggregates for both candidates. Cubic improves 3-D RMS on all ten seeds in every one of the eight stationary scenarios. Fitted does the same.

By spectral family, cubic changes 3-D RMS by -4.84% on JONSWAP and -4.25% on PM-Stokes. Fitted changes it by -5.79% on JONSWAP and -4.73% on PM-Stokes.

### Per-scenario stationary 3-D RMS change

| sea | cubic p=3 | fitted p=2.9052 | seeds improved |
|---|---:|---:|---:|
| JONSWAP Hs=0.27 m | -11.23% | -11.16% | 10/10 both |
| JONSWAP Hs=1.50 m | -9.28% | -9.33% | 10/10 both |
| JONSWAP Hs=4.00 m | -5.92% | -6.40% | 10/10 both |
| JONSWAP Hs=8.50 m | -3.38% | -4.71% | 10/10 both |
| PM-Stokes Hs=0.27 m | -10.76% | -10.71% | 10/10 both |
| PM-Stokes Hs=1.50 m | -7.86% | -7.90% | 10/10 both |
| PM-Stokes Hs=4.00 m | -5.41% | -5.64% | 10/10 both |
| PM-Stokes Hs=8.50 m | -2.91% | -3.60% | 10/10 both |

The corresponding cubic Z-error changes in percentage points of Hs are -0.0746, +0.0206, +0.0882, +0.0292 for the four JONSWAP seas and -0.0382, +0.1352, +0.1040, +0.0715 for PM-Stokes. Thus the stationary 3-D improvement is driven by the horizontal channels and is not a uniform all-axis improvement.

## Cubic versus exact fitted exponent

Comparing fitted p=2.9052 directly against cubic p=3 on the eight stationary seas:

- fitted lowers 3-D RMS a further 0.764% [0.538%, 1.015%] relative to cubic, with the same sign on all ten seed aggregates;
- fitted lowers Y RMS another 1.72%;
- fitted lowers X RMS another 0.28%;
- fitted raises Z RMS by 0.509% and raises Z error by 0.01025 percentage points of Hs relative to cubic, with the same sign on all ten seed aggregates.

So the exact exponent buys a small additional horizontal/3-D improvement at a measurable additional vertical cost.

## Non-stationary transition

For the complete transition scoring window, cubic changes 3-D RMS by -0.218% (6/10 seed triplets improve) while increasing Z RMS by +4.26% (0/10 improve). Fitted changes 3-D RMS by -0.272% (6/10 improve) while increasing Z RMS by +5.54% (0/10 improve).

Segmented transition results:

| segment | cubic Z | cubic 3-D | fitted Z | fitted 3-D |
|---|---:|---:|---:|---:|
| start | +0.86% | -9.82% | +0.88% | -9.86% |
| blend | +8.52% | +3.62% | +10.08% | +4.28% |
| end | +2.13% | -1.53% | +3.30% | -1.89% |

For cubic, the blend-segment Z degradation occurs on all 10 seeds and the blend 3-D degradation occurs on 9/10. For fitted, the corresponding blend results are also 10/10 Z worse and 9/10 3-D worse. At the end segment both candidates improve 3-D on all 10 seeds but worsen Z on all 10.

## Interpretation

The complete validation supports the physical-period result that a near-cubic effective law is much better for the stationary full-3D objective than the deployed tau^2.5 law. However, the nominal full-3D normalization found on the stationary calibration sea is too aggressive to be called a universal production optimum: it systematically sacrifices vertical accuracy and introduces a clear transient penalty during the transition blend.

If the sole objective is stationary full-3D displacement RMS, fitted p=2.9052 is best of these three tested laws. If a simple law is desired, cubic p=3 captures most of that benefit with less vertical and transition damage. Neither candidate dominates the deployed law on all required metrics, so neither should replace the deployed default yet without a normalization / transient-scheduling follow-up experiment.
