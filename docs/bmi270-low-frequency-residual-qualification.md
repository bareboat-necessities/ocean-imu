# AtomS3R/BMI270 low-frequency residual qualification

This protocol qualifies the assembled device for the OU-III finite-error proof.
It does not change shipping filtering, calibration, tuning, or estimator behavior.

## Bound quantity

Use the actual shipping BMI270 configuration, calibration/temperature
compensation, mounting, power path, timestamps and units.  Form residuals
against an independent reference:

    n_a = a_measured_calibrated - a_reference_specific_force
    n_g = gyro_measured_calibrated - gyro_reference

Apply the fixed offline certification operator L_X to each residual vector.
The qualification values include reference uncertainty:

    eps_a_LF = max ||L_X n_a|| + u_a_ref
    eps_g_LF = max ||L_X n_g|| + u_g_ref

No subtraction of reference uncertainty is allowed.

## L_X

The qualification operator is a zero-phase offline low-pass.  Its required
response is:

- passband: 0--0.15 Hz, gain >= 0.99;
- nominal transition: 0.20--0.30 Hz;
- stopband begins at 0.30 Hz;
- zero phase (offline forward/backward or symmetric FIR).

The qualification artifact MUST export the actual coefficients, sampled
frequency response, implementation version and hash.  The proof binds those
artifacts rather than the phrase "0.25-Hz filter".

## Required assembled-device captures

### Stationary thermal/orientation

Capture approximately +X,-X,+Y,-Y,+Z,-Z orientations.  Angular-rate reference
is zero; reference specific force is local gravity in the surveyed
orientation.  Repeat across the declared deployment temperature range,
multiple power cycles and representative mounting/PCB stress.  Records must
be long enough to expose behavior below 0.02 Hz; use hours, not minutes.

### Slow single-axis dynamics

Using an independent encoder/rate-table reference, exercise roll and pitch at

    .02, .05, .08, .10, .15, .20 Hz

and amplitudes spanning approximately 1--10 degrees.  Use many cycles per
condition.  Compute reference gravity/specific force and angular rate at the
IMU location, including the configured lever arm where applicable.

### Combined roll/pitch

Repeat representative two-axis trajectories in the same frequency/amplitude
range.  Qualification is on vector norm, not independent axis maxima.

### Temperature/repeatability

Repeat dynamic conditions over temperature and after multiple power cycles.
The theorem qualifies the assembled device population/envelope, not a bare
BMI270 typical specification.

## Capture CSV contract

One row per shipping IMU sample:

    t_s,temp_C,
    ax_mps2,ay_mps2,az_mps2,
    gx_rad_s,gy_rad_s,gz_rad_s,
    ax_ref_mps2,ay_ref_mps2,az_ref_mps2,
    gx_ref_rad_s,gy_ref_rad_s,gz_ref_rad_s

Metadata stored beside each capture must include device/PCB ID, firmware and
calibration hashes, mounting, power source, BMI270 register/configuration
fingerprint, reference instrument and calibration, orientation/motion command,
temperature condition and uncertainty bounds u_a_ref/u_g_ref.

## 1 degree / 60 s proof gate

For T_X=60 s and theta_X=1 degree, after charging the existing physical
bias-rate ambiguity, the strict residual qualification is

    2 asin(eps_a_LF / 9.80665) + 60 eps_g_LF
      < 0.0113349952420757 rad.

This is a JOINT tradeoff.  Do not independently require both axis intercepts.

Axis intercepts:

    eps_a_LF < .0555788680 m/s2       (negligible LF gyro residual)
    eps_g_LF < 1.88916587e-4 rad/s    (negligible LF accel residual)

A convenient initial 50/50 laboratory target is approximately

    eps_a_LF < .02779 m/s2
    eps_g_LF < 9.45e-5 rad/s

with additional qualification margin.

## Publication artifact

Publish the worst assembled-device pair (eps_a_LF,eps_g_LF), joint angular
charge and strict margin, per-capture maxima, reference uncertainty, device and
condition coverage, L_X coefficients/response/hash, and raw-data hashes.

Passing finite captures is evidence for the declared deployment qualification;
it is not by itself a mathematical all-history theorem.  The theorem must state
the resulting deterministic low-frequency envelope explicitly.
