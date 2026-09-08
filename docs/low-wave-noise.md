# Low-wave acceleration weighting

The OU-II and OU-III simulators configure a known 28 ft vessel response. The generic filters keep this optional measurement weighting disabled. It is independent of the direction equalizer switch and does not change the OU time constant, stationary covariance, pseudo-measurement coefficients, sensor noise injection, or acceptance thresholds.

For each horizontal sensor axis, the preceding measurement-only wave-band standard deviation and frequency predict a conservative horizontal signal level:

\[
\hat\sigma_i=\sigma_z\frac{\sqrt{(1-(fT_h)^2)^2+(2\zeta_h fT_h)^2}}{1+(2\pi f\tau_i)^2},
\qquad \rho_i=\hat\sigma_i/\sigma_{a,i}.
\]

The unknown incident-direction projection is bounded by one. Common footprint/depth attenuation is omitted conservatively. This scalar approximation is a regularization policy for the configured stationary hull, not an inferred sensor-noise density or a general directional RAO inversion.

The accelerometer measurement standard-deviation multiplier is four for SNR at most two and exactly one for SNR at least four. Between those points, its variance multiplier follows a cubic smoothstep from 16 to one. Vertical measurement weighting is unchanged. Vibration inflation is added in variance after this weighting. Before the adaptive filter becomes live, the existing startup covariance applies. The current sample cannot choose its own weighting: frequency and wave-band scale are read before that sample updates the tuner.

`SF_LOW_WAVE_RACC_MAX_SCALE=1` disables the simulation policy. `SF_LOW_WAVE_RACC_SNR` sets the lower transition endpoint (default two; the upper endpoint is twice that value). `SF_SIGMA_A_X_SCALE`, `_Y_SCALE`, and `_Z_SCALE` provide explicit static measurement-weight ablations. These controls never alter injected noise.

The policy improves acceleration noise near the response floor, with small attitude tradeoffs. It does not resolve the static roll/accelerometer-bias ambiguity, guarantee wave-direction availability, or calibrate magnetic distortion. In low waves, a longer OU time constant alone does not imply smoother acceleration: process driving, covariance synchronization, and measurement/pseudo-measurement weights also determine the estimator bandwidth. The retained OU-III R_S horizontal standard-deviation factors remain 0.72; its stationary horizontal/vertical prior factor remains one. OU-II retains its prior factor of 1.5 and MSE ratio of 0.3.

The paired study, rejected alternatives, fresh sensor draws, heading/incident-angle checks, and signal-band response measurements are recorded in `reports/results/low_wave_noise/`. Angular errors are conditional on an available axis; availability and unresolved travel fractions must be read alongside them. OU-III applied-parameter telemetry reads the active MEKF matrices, including floors and channel freezes, rather than the next staged tuner command.
