# Physical-period r_S adaptation optimum

This result is from the successful refinement of the dedicated full-MEKF experiment on draft PR #324, GitHub Actions run 31916274242 attempt 2 (job 95088661303), commit c869be0dbb6a12502f42bc51a5daf8c3782edb69.

## Question

Measure the adaptation exponent with respect to *physical sea period*, rather than perturbing the OU filter parameter on a fixed physical sea.

The source is the JONSWAP Hs=1.5 m, L=50.710 m stationary record. Five physical time stretches s = {0.85, 0.925, 1.0, 1.075, 1.15} are constructed. Translational acceleration RMS is held fixed exactly (to numerical precision) while displacement and velocity obey p_s(t)=c s^2 p(t/s), v_s(t)=c s v(t/s), a_s(t)=c a(t/s). Attitude is time-stretched at fixed angular standard deviation and the body IMU is rebuilt, so attitude/gravity leakage remains present.

At each stretch, tau_filter = tau_0 s, sigma_aw is fixed, and the effective r_S standard deviation actually entering Kalman3D_Wave_OU_III is swept independently. The wrapper's cadence information-rate factor is explicitly inverted when constructing the fixed-tuning command, so it cannot predetermine the measured exponent.

Three matched wave/IMU/initialization seed triplets are used. Each replay is 600 s and the trailing 300 s is scored. The refined grid contains nine effective-r_S factors {0.25, 0.35, 0.45, 0.55, 0.65, 0.78, 0.9, 1.05, 1.2} around the deployed tau^2.5 centre, for 135 full-MEKF replays. At every stretch the optimum is interior to the refined grid.

## Result

For vertical displacement MSE,

- p_tau = **3.0053**
- seed-triplet bootstrap 95% descriptive interval: **2.6517 to 3.2844**

For full 3-D displacement MSE,

- p_tau = **2.9052**
- seed-triplet bootstrap 95% descriptive interval: **2.8346 to 3.0276**

The bootstrap interval is descriptive because there are only three independent seed triplets; the main evidence is the smooth interior optimum at all five physical-period points.

### Vertical-MSE optima

| physical stretch | effective r_S optimum [m s] | deployed tau^2.5 centre [m s] | optimum/deployed |
|---:|---:|---:|---:|
| 0.850 | 1.025892 | 1.241597 | 0.8263 |
| 0.925 | 1.325828 | 1.533867 | 0.8644 |
| 1.000 | 1.578124 | 1.863947 | 0.8467 |
| 1.075 | 2.045023 | 2.233339 | 0.9157 |
| 1.150 | 2.575293 | 2.643492 | 0.9742 |

### Full-3-D-MSE optima

| physical stretch | effective r_S optimum [m s] | deployed tau^2.5 centre [m s] | optimum/deployed |
|---:|---:|---:|---:|
| 0.850 | 0.741031 | 1.241597 | 0.5968 |
| 0.925 | 0.941594 | 1.533867 | 0.6139 |
| 1.000 | 1.160579 | 1.863947 | 0.6226 |
| 1.075 | 1.448576 | 2.233339 | 0.6486 |
| 1.150 | 1.790962 | 2.643492 | 0.6775 |

## Local fitted laws

Using the nominal calibration point tau_0 = 2.17904091 s and sigma_0 = 0.724445343 m/s^2, and combining this period experiment with the separately measured amplitude exponent p_sigma = 1, the fitted effective-input laws are

- vertical objective: r_S,eff ~= **0.20968 sigma_aw tau^3.0053**
- full 3-D objective: r_S,eff ~= **0.16671 sigma_aw tau^2.9052**

The currently deployed wrapper instead produces r_S,eff ~= 0.36708 sigma_aw tau^2.5 at the filter input.

If the current cadence normalization r_S,eff = r_S,base sqrt(T_S0/T_S), with T_S proportional to tau, is retained, then the equivalent *base-command* laws are

- vertical objective: r_S,base ~= **0.19992 sigma_aw tau^3.5053**
- full 3-D objective: r_S,base ~= **0.15895 sigma_aw tau^3.4052**

Equivalently, removing the cadence information-rate renormalization and applying the effective law directly would avoid the extra +1/2 exponent in the base command.

## Interpretation

This local full-MEKF experiment rejects tau^2.5 as the optimum period dependence around the nominal JONSWAP sea. The vertical optimum is essentially cubic in physical period; the full-3-D optimum is slightly shallower but still close to cubic. The result is local to a self-similar period family around the nominal sea and should be confirmed on additional spectral families before being promoted to a universal adaptation law.
