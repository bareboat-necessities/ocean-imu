In simple terms:

- Common Quaternion Multiplicative EKF (Extended Kalman Filter) model
model with additional states for 3D kinematics
- Acceleration measurement regularized with OU stochastic autocorrelation process
- Whole thing in contious time is LTE stochastic
- Over dt step is solved via matrix exponential
- With small dt matrices can become numericsllly ill conditioned
- So Joseph form and small series are used
- Double integral drift correction measures
- Adaptation of Kalman parameters required to match sea all sea states
- Adapration is using acceleratetion freq tracking for time scale and acceleration variance direct estimatiion for amplitude scale
- Making sure timestamps of samples from IMU are correct using FIFO timestamps
- 
