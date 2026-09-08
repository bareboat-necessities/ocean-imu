# Wave direction with a known vessel RAO

The OU-II and OU-III direction branches use the engine-guarded accelerometer, level it using the attitude estimate, and optionally match the three translational RAOs before estimating the propagation axis and travel sense. The wave-period tuner reads the pre-equalized signal. Equalization does not feed back into attitude, displacement, or adaptive OU regularization. TFG currently has no direction estimator.

## Complex response matching

Amplitude correction alone is insufficient: unequal surge/sway phase changes horizontal covariance, and heave phase determines the vertical–horizontal quadrature used for travel sense. For the stationary 28 ft preset, using the causal Laplace convention,

\[
Q_x(s)=(1+0.7s)^{-2},\qquad Q_y(s)=(1+1.0s)^{-2},\qquad
G_h(s)=\left(1+2\zeta_hs/\omega_h+s^2/\omega_h^2\right)^{-1},
\]

where \(\omega_h=2\pi/2.4\) and \(\zeta_h=0.45\). Horizontal motion also contains the common footprint/depth response and the incident directional components. Match the channels to \(Q_t=(1+\tau_t s)^{-2}\):

\[
E_x=Q_t/Q_x,\qquad E_y=Q_t/Q_y,\qquad E_z=Q_t/G_h.
\]

For this preset \(\tau_t=1\) second. Each compensator is stable, proper, has unity DC gain, and has magnitude at most one. This corrects relative amplitude and phase while retaining common low-pass attenuation. An inverse to unity would amplify high-frequency measurement noise. Bilinear biquads implement these responses; an independent complex-harmonic unit test covers five frequencies and five directions, including above the heave resonance.

The generator uses `Re[H exp(i phi - i omega t)]`; consequently its printed complex transfer functions are conjugates of the positive-frequency causal convention above. Its `cos(k·x - omega t + phi)` phase kernel makes filename azimuth the propagation-to direction. Axial scoring remains modulo 180 degrees; directed scoring uses that propagation-to reference modulo 360 degrees.

## Coordinates and scope

The generator's physical axes are X forward, Y port, Z up. The simulation's existing up-to-NED conversion swaps the horizontal axes. At the fixed zero heading in these records, the direction stage therefore receives sway in its first horizontal channel and surge in its second; its configured time constants are 1.0 and 0.7 seconds respectively. The existing direction output mapping converts back to generator azimuth. A deployed hull must configure the response in its actual sensor/heading basis rather than copy these simulation channel assignments blindly.

The generic equalizer is disabled by default. Both OU simulators explicitly enable `vessel-rao-28ft`; `W3D_DIRECTION_RAO=off` provides a paired ablation. The profile assumes the known stationary hull. Forward speed, heading-dependent response, nonzero sensor lever arm, and another hull require the appropriate response model. A general directional RAO may require fitting a directional spectrum through the complex transfer matrix instead of this separable three-channel equalizer.

Common footprint/depth attenuation is retained. This is not incident-wave-height reconstruction, nor can it create information at response nulls or resolve every broad directional sea. Existing confidence and acceptance thresholds remain unchanged. Weak-motion cases can still withhold travel direction, and attenuation can reduce their availability.

## Validation

`tools/direction_rao_ablation.py` runs all eight records with OU-II and OU-III, engine off and nominal 2400 RPM (0.6 m/s² reference vibration, 80 Hz sensor bandwidth), and three arms: unconditioned, guard only, and guard plus RAO. The 96 full-record replays use the final 900 seconds and retain all quality failures. Results and exact input/binary hashes are under `reports/results/direction_rao_ablation/`.

Guard-off disables both conditioning and vibration covariance inflation. Quiet guard-on/off parity checks detector transparency. Guard-only versus guard-plus-RAO pairs isolate direction equalization; their upstream attitude and displacement measurements must be identical. Report unresolved fractions alongside angular errors so improved accuracy cannot hide reduced availability.
