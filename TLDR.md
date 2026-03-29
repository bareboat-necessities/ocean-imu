In simple terms:

- Common Quaternion Multiplicative EKF (Extended Kalman Filter) model, with additional states for 3D kinematics in ocean waves. 
- Acceleration measurements are regularized with an Ornstein–Uhlenbeck (OU) stochastic autocorrelation process.
- The whole continuous-time system is linear time-varying (LTV) stochastic.
- Each `dt` step is solved via matrix exponential.
- With very small `dt`, matrices can become numerically ill-conditioned (a.k.a. math gets seasick).
- So Joseph form and small-series approximations are used for numerical stability.
- Includes double-integration drift-correction pseudo-measurements.
- Kalman parameters must be adapted to match actual sea states.
- Adaptation uses acceleration-frequency tracking for time-scale tuning, and direct acceleration-variance estimation for amplitude-scale tuning.
- IMU sample timestamps are using accurate FIFO timestamps. Correct sample timestamps are very important due to samples being double integrated over time steps.
- Performs initial learning of gravity relation to magnetic field needed because sensor might be turned on when already in waves motion
