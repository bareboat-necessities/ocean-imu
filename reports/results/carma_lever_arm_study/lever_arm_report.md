# CARMA IMU lever-arm installation study

The wave records define motion at the vessel centre of gravity (CG).  An IMU
rigidly mounted at body offset $r$ additionally measures
`alpha x r + omega x (omega x r)`; the simulator applies that term to the
accelerometer truth before sensor corruption.  Nothing else changes between
arms: noise realization, CARMA configuration, startup, adaptation and scoring
are those of the gated `carma3-sim` (A41) and `tfg_carma3-sim` (T41) replays.

| Arm | Correction |
| --- | --- |
| baseline | none needed: IMU at the CG |
| unmodeled | none |
| exact | simulator removes the installed term from truth kinematics (oracle) |
| wrapper | `setImuLeverArmBody`, the deployed correction: measured rate, 15 Hz two-pole low-pass, second-order difference, at the wrapper input |
| mekf | ablation: the same model inside the inner filter's accelerometer update |
| mekf_ou3 | ablation: OU-III's model (one-pole first difference, alpha bias Jacobian) inside the inner filter |

Offsets are in the sketches' body frame: forward, starboard, and above the CG.
Scoring uses the trailing **900 s** of each 1200 s record; every
entry is the RMS pooled over the eight seas; ratios are to the filter's own
CG baseline; the worst-sea columns divide each sea by its own baseline and
keep the largest.  The inner-filter arms report no residual: the simulator cannot
see a correction made inside the estimator.  The removed fraction of the
uncorrected 3-D excess is shown where that excess exceeds 5 % of the baseline.

Source commit: `1201f2c15274ee95fc230c77a280c2d8ec8ee47f`.

## CARMA3

CG baseline: 3-D displacement 0.2371 m, max roll/pitch 0.2056 deg, yaw 0.7137 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.001x | 1.053x | 0.963x | 0.996x | 0.996x | 0.0161 | 0.0161 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0161 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 0.998x | 1.007x | 1.010x | 0.0161 | n/a | n/a |
| mekf | x-forward | 10 | 1.001x | 1.003x | 1.002x | 1.002x | 1.008x | 0.0161 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.003x | 1.004x | 1.105x | 1.111x | 0.998x | 0.0161 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.130x | 0.931x | 0.996x | 0.992x | 0.0322 | 0.0322 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0322 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.003x | 0.997x | 1.013x | 1.020x | 0.0322 | n/a | n/a |
| mekf | x-forward | 20 | 1.004x | 1.012x | 1.007x | 1.025x | 1.018x | 0.0322 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.007x | 1.013x | 1.225x | 1.234x | 0.997x | 0.0322 | n/a | n/a |
| unmodeled | x-forward | 30 | 0.972x | 1.226x | 0.900x | 0.992x | 0.985x | 0.0483 | 0.0483 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0483 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.005x | 0.996x | 1.019x | 1.028x | 0.0483 | n/a | n/a |
| mekf | x-forward | 30 | 0.970x | 1.011x | 1.011x | 1.009x | 1.023x | 0.0483 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 0.976x | 1.014x | 1.353x | 1.361x | 1.010x | 0.0483 | n/a | n/a |
| unmodeled | x-forward | 50 | 0.970x | 1.456x | 0.925x | 0.990x | 0.987x | 0.0805 | 0.0805 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0805 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.009x | 0.994x | 1.022x | 1.051x | 0.0805 | n/a | n/a |
| mekf | x-forward | 50 | 0.962x | 1.009x | 1.017x | 1.026x | 1.050x | 0.0805 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 0.974x | 1.022x | 1.625x | 1.721x | 1.024x | 0.0805 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.011x | 2.103x | 1.026x | 1.112x | 1.006x | 0.1609 | 0.1609 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.1609 | 0.0000 | n/a |
| wrapper | x-forward | 100 | 1.000x | 1.017x | 0.989x | 1.035x | 1.137x | 0.1609 | n/a | n/a |
| mekf | x-forward | 100 | 0.967x | 1.025x | 1.029x | 1.046x | 1.140x | 0.1609 | n/a | n/a |
| mekf_ou3 | x-forward | 100 | 1.006x | 1.077x | 2.687x | 2.950x | 2.139x | 0.1609 | n/a | n/a |
| unmodeled | x-forward | 200 | 1.161x | 3.363x | 1.261x | 1.557x | 1.107x | 0.3218 | 0.3218 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3218 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.007x | 1.062x | 0.982x | 1.047x | 1.165x | 0.3218 | n/a | 95% |
| mekf | x-forward | 200 | 0.980x | 1.034x | 1.019x | 1.260x | 1.477x | 0.3218 | n/a | 113% |
| mekf_ou3 | x-forward | 200 | 1.232x | 1.506x | 5.649x | 6.317x | 8.143x | 0.3218 | n/a | -44% |
| unmodeled | y-starboard | 10 | 1.004x | 1.021x | 0.989x | 1.023x | 1.001x | 0.0087 | 0.0087 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0087 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.002x | 1.005x | 0.991x | 0.0087 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.011x | 1.002x | 1.006x | 0.990x | 0.0087 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 1.009x | 1.016x | 1.015x | 0.954x | 0.0087 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.008x | 1.043x | 0.975x | 1.061x | 1.007x | 0.0173 | 0.0173 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0173 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.005x | 1.009x | 0.984x | 0.0173 | n/a | n/a |
| mekf | y-starboard | 20 | 0.999x | 1.016x | 1.002x | 1.013x | 0.987x | 0.0173 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.996x | 1.013x | 1.029x | 1.033x | 0.944x | 0.0173 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.014x | 1.061x | 0.964x | 1.062x | 1.009x | 0.0260 | 0.0260 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0260 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.006x | 1.015x | 0.982x | 0.0260 | n/a | n/a |
| mekf | y-starboard | 30 | 0.999x | 1.012x | 1.003x | 1.019x | 0.981x | 0.0260 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.995x | 1.008x | 1.046x | 1.054x | 0.963x | 0.0260 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.030x | 1.160x | 0.998x | 1.091x | 1.015x | 0.0434 | 0.0434 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0434 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.000x | 1.001x | 1.009x | 1.022x | 0.983x | 0.0434 | n/a | n/a |
| mekf | y-starboard | 50 | 0.999x | 1.007x | 1.005x | 1.035x | 0.973x | 0.0434 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.995x | 1.004x | 1.085x | 1.107x | 1.102x | 0.0434 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.090x | 1.591x | 1.132x | 1.193x | 1.054x | 0.0867 | 0.0867 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0867 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 0.997x | 1.001x | 1.007x | 1.030x | 1.014x | 0.0867 | n/a | 103% |
| mekf | y-starboard | 100 | 1.002x | 1.012x | 0.997x | 1.035x | 0.995x | 0.0867 | n/a | 98% |
| mekf_ou3 | y-starboard | 100 | 1.011x | 1.025x | 1.864x | 2.290x | 2.014x | 0.0867 | n/a | 87% |
| unmodeled | y-starboard | 200 | 1.233x | 2.725x | 1.321x | 1.400x | 1.124x | 0.1734 | 0.1734 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1734 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.002x | 1.025x | 0.983x | 0.997x | 1.022x | 0.1734 | n/a | 99% |
| mekf | y-starboard | 200 | 0.967x | 1.021x | 0.966x | 1.009x | 1.244x | 0.1734 | n/a | 114% |
| mekf_ou3 | y-starboard | 200 | 1.191x | 1.321x | 6.681x | 7.654x | 7.185x | 0.1734 | n/a | 18% |
| unmodeled | z-up | 10 | 1.009x | 1.035x | 0.999x | 1.003x | 0.999x | 0.0181 | 0.0181 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0181 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.006x | 0.990x | 0.0181 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.000x | 1.003x | 1.008x | 0.990x | 0.0181 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 0.999x | 1.003x | 1.006x | 0.993x | 0.0181 | n/a | n/a |
| unmodeled | z-up | 20 | 1.021x | 1.098x | 0.998x | 1.017x | 0.998x | 0.0362 | 0.0362 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0362 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.001x | 1.005x | 1.012x | 0.985x | 0.0362 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.004x | 1.006x | 1.017x | 0.981x | 0.0362 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.999x | 1.002x | 1.005x | 1.012x | 0.985x | 0.0362 | n/a | n/a |
| unmodeled | z-up | 30 | 1.034x | 1.184x | 0.997x | 1.021x | 0.997x | 0.0543 | 0.0543 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0543 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.000x | 1.001x | 1.008x | 1.020x | 0.981x | 0.0543 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.006x | 1.010x | 1.025x | 0.976x | 0.0543 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.002x | 1.008x | 1.019x | 0.981x | 0.0543 | n/a | n/a |
| unmodeled | z-up | 50 | 1.066x | 1.409x | 0.996x | 1.023x | 0.996x | 0.0904 | 0.0904 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0904 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.000x | 1.002x | 1.014x | 1.032x | 0.988x | 0.0904 | n/a | 99% |
| mekf | z-up | 50 | 1.001x | 1.005x | 1.017x | 1.045x | 0.984x | 0.0904 | n/a | 98% |
| mekf_ou3 | z-up | 50 | 0.995x | 1.000x | 1.015x | 1.034x | 0.973x | 0.0904 | n/a | 107% |
| unmodeled | z-up | 100 | 1.179x | 2.130x | 0.991x | 1.030x | 0.989x | 0.1808 | 0.1808 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1808 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.004x | 1.021x | 1.044x | 1.148x | 0.1808 | n/a | 99% |
| mekf | z-up | 100 | 1.003x | 1.006x | 1.028x | 1.083x | 1.158x | 0.1808 | n/a | 99% |
| mekf_ou3 | z-up | 100 | 0.991x | 1.002x | 1.030x | 1.072x | 0.996x | 0.1808 | n/a | 105% |
| unmodeled | z-up | 200 | 1.500x | 3.791x | 0.986x | 1.068x | 0.991x | 0.3617 | 0.3617 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.3617 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.008x | 1.028x | 1.012x | 1.192x | 1.604x | 0.3617 | n/a | 98% |
| mekf | z-up | 200 | 1.020x | 1.089x | 1.234x | 2.780x | 3.615x | 0.3617 | n/a | 96% |
| mekf_ou3 | z-up | 200 | 0.987x | 1.023x | 1.056x | 1.122x | 1.327x | 0.3617 | n/a | 103% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2372 m, max roll/pitch 0.2056 deg, yaw 0.7128 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.000x | 1.053x | 0.963x | 0.996x | 0.994x | 0.0161 | 0.0161 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0161 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 0.999x | 1.003x | 1.003x | 0.0161 | n/a | n/a |
| mekf | x-forward | 10 | 1.001x | 1.003x | 1.003x | 1.002x | 1.002x | 0.0161 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.003x | 1.004x | 1.108x | 1.114x | 0.994x | 0.0161 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.130x | 0.931x | 0.994x | 0.990x | 0.0322 | 0.0322 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0322 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.003x | 0.998x | 1.007x | 1.008x | 0.0322 | n/a | n/a |
| mekf | x-forward | 20 | 1.004x | 1.012x | 1.009x | 1.018x | 1.006x | 0.0322 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.007x | 1.013x | 1.230x | 1.238x | 0.992x | 0.0322 | n/a | n/a |
| unmodeled | x-forward | 30 | 0.972x | 1.226x | 0.900x | 0.990x | 0.983x | 0.0483 | 0.0483 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0483 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.005x | 0.998x | 1.011x | 1.011x | 0.0483 | n/a | n/a |
| mekf | x-forward | 30 | 0.970x | 1.011x | 1.014x | 1.010x | 1.005x | 0.0483 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 0.975x | 1.013x | 1.361x | 1.369x | 1.004x | 0.0483 | n/a | n/a |
| unmodeled | x-forward | 50 | 0.970x | 1.456x | 0.926x | 0.987x | 0.981x | 0.0805 | 0.0805 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0805 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.009x | 0.997x | 1.012x | 1.021x | 0.0805 | n/a | n/a |
| mekf | x-forward | 50 | 0.962x | 1.009x | 1.020x | 1.019x | 1.015x | 0.0805 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 0.973x | 1.022x | 1.632x | 1.726x | 0.998x | 0.0805 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.011x | 2.103x | 1.026x | 1.120x | 0.999x | 0.1609 | 0.1609 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.1609 | 0.0000 | n/a |
| wrapper | x-forward | 100 | 1.000x | 1.017x | 0.993x | 1.025x | 1.055x | 0.1609 | n/a | n/a |
| mekf | x-forward | 100 | 0.967x | 1.025x | 1.034x | 1.029x | 1.047x | 0.1609 | n/a | n/a |
| mekf_ou3 | x-forward | 100 | 1.001x | 1.077x | 2.390x | 3.089x | 1.627x | 0.1609 | n/a | n/a |
| unmodeled | x-forward | 200 | 1.161x | 3.363x | 1.266x | 1.570x | 1.102x | 0.3218 | 0.3218 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.3218 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.006x | 1.062x | 0.985x | 1.058x | 1.060x | 0.3218 | n/a | 96% |
| mekf | x-forward | 200 | 0.976x | 1.034x | 1.028x | 1.136x | 1.168x | 0.3218 | n/a | 115% |
| mekf_ou3 | x-forward | 200 | 1.111x | 1.327x | 3.764x | 8.100x | 11.874x | 0.3218 | n/a | 31% |
| unmodeled | y-starboard | 10 | 1.004x | 1.021x | 0.989x | 1.023x | 1.000x | 0.0087 | 0.0087 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0087 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.001x | 1.002x | 0.994x | 0.0087 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.011x | 1.001x | 1.004x | 0.993x | 0.0087 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 1.009x | 1.015x | 1.015x | 0.958x | 0.0087 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.008x | 1.043x | 0.975x | 1.062x | 1.007x | 0.0173 | 0.0173 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0173 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.002x | 1.004x | 0.990x | 0.0173 | n/a | n/a |
| mekf | y-starboard | 20 | 0.999x | 1.016x | 0.999x | 1.008x | 0.993x | 0.0173 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.996x | 1.013x | 1.027x | 1.030x | 0.952x | 0.0173 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.014x | 1.061x | 0.964x | 1.062x | 1.007x | 0.0260 | 0.0260 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0260 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.003x | 1.006x | 0.989x | 0.0260 | n/a | n/a |
| mekf | y-starboard | 30 | 0.999x | 1.012x | 1.000x | 1.012x | 0.985x | 0.0260 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.995x | 1.008x | 1.042x | 1.048x | 0.970x | 0.0260 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.030x | 1.160x | 1.003x | 1.091x | 1.016x | 0.0434 | 0.0434 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0434 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.000x | 1.000x | 1.004x | 1.009x | 0.986x | 0.0434 | n/a | n/a |
| mekf | y-starboard | 50 | 0.999x | 1.007x | 1.000x | 1.025x | 0.977x | 0.0434 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.995x | 1.003x | 1.077x | 1.083x | 1.113x | 0.0434 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.090x | 1.591x | 1.139x | 1.188x | 1.056x | 0.0867 | 0.0867 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0867 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 0.997x | 1.000x | 1.002x | 1.010x | 0.989x | 0.0867 | n/a | 103% |
| mekf | y-starboard | 100 | 1.002x | 1.012x | 0.994x | 1.024x | 0.978x | 0.0867 | n/a | 98% |
| mekf_ou3 | y-starboard | 100 | 1.009x | 1.036x | 1.723x | 2.922x | 1.936x | 0.0867 | n/a | 90% |
| unmodeled | y-starboard | 200 | 1.232x | 2.725x | 1.328x | 1.418x | 1.125x | 0.1734 | 0.1734 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.1734 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.001x | 1.025x | 0.988x | 1.001x | 0.993x | 0.1734 | n/a | 99% |
| mekf | y-starboard | 200 | 0.965x | 1.017x | 0.970x | 0.989x | 1.019x | 0.1734 | n/a | 115% |
| mekf_ou3 | y-starboard | 200 | 1.124x | 3.612x | 4.187x | 9.337x | 40.147x | 0.1734 | n/a | 47% |
| unmodeled | z-up | 10 | 1.009x | 1.035x | 1.000x | 1.003x | 0.998x | 0.0181 | 0.0181 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0181 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.001x | 1.003x | 0.993x | 0.0181 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.005x | 0.992x | 0.0181 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 0.999x | 1.002x | 1.004x | 0.994x | 0.0181 | n/a | n/a |
| unmodeled | z-up | 20 | 1.021x | 1.098x | 0.999x | 1.017x | 0.999x | 0.0362 | 0.0362 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0362 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.001x | 1.002x | 1.004x | 0.988x | 0.0362 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.004x | 1.003x | 1.011x | 0.988x | 0.0362 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.999x | 1.002x | 1.004x | 1.009x | 0.991x | 0.0362 | n/a | n/a |
| unmodeled | z-up | 30 | 1.034x | 1.185x | 1.000x | 1.019x | 0.995x | 0.0543 | 0.0543 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0543 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.000x | 1.001x | 1.003x | 1.008x | 0.984x | 0.0543 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.006x | 1.005x | 1.015x | 0.979x | 0.0543 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.002x | 1.006x | 1.013x | 0.983x | 0.0543 | n/a | n/a |
| unmodeled | z-up | 50 | 1.066x | 1.409x | 1.000x | 1.020x | 0.994x | 0.0904 | 0.0904 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.0904 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.000x | 1.001x | 1.006x | 1.013x | 0.976x | 0.0904 | n/a | 100% |
| mekf | z-up | 50 | 1.001x | 1.005x | 1.009x | 1.027x | 0.970x | 0.0904 | n/a | 98% |
| mekf_ou3 | z-up | 50 | 0.995x | 1.000x | 1.010x | 1.023x | 0.975x | 0.0904 | n/a | 107% |
| unmodeled | z-up | 100 | 1.178x | 2.131x | 0.999x | 1.030x | 0.985x | 0.1808 | 0.1808 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1808 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.001x | 1.003x | 1.010x | 1.024x | 1.010x | 0.1808 | n/a | 99% |
| mekf | z-up | 100 | 1.002x | 1.005x | 1.017x | 1.057x | 1.000x | 0.1808 | n/a | 99% |
| mekf_ou3 | z-up | 100 | 0.991x | 1.002x | 1.021x | 1.049x | 0.956x | 0.1808 | n/a | 105% |
| unmodeled | z-up | 200 | 1.499x | 3.791x | 1.001x | 1.069x | 0.980x | 0.3617 | 0.3617 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 0.999x | 0.3617 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.006x | 1.028x | 1.005x | 1.042x | 1.146x | 0.3617 | n/a | 99% |
| mekf | z-up | 200 | 1.008x | 1.028x | 1.028x | 1.084x | 1.556x | 0.3617 | n/a | 98% |
| mekf_ou3 | z-up | 200 | 0.985x | 1.023x | 1.043x | 1.105x | 1.045x | 0.3617 | n/a | 103% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.023x | 1.007x | 0.997x |
| CARMA3 | 2 | 1.005x | 1.004x | 0.997x |
| CARMA3 | 5 | 1.001x | 1.005x | 0.988x |
| CARMA3 | 10 | 1.001x | 1.014x | 0.994x |
| CARMA3 | 15 | 1.002x | 1.021x | 1.148x |
| CARMA3 | 25 | 1.007x | 1.007x | 1.406x |
| CARMA3 | 50 | 1.025x | 0.988x | 0.970x |
| CARMA3 | 100 | 1.041x | 0.988x | 0.957x |
| TFG-CARMA3 | 1 | 1.023x | 1.007x | 0.997x |
| TFG-CARMA3 | 2 | 1.005x | 1.003x | 0.997x |
| TFG-CARMA3 | 5 | 1.001x | 1.002x | 0.992x |
| TFG-CARMA3 | 10 | 1.001x | 1.006x | 0.978x |
| TFG-CARMA3 | 15 | 1.001x | 1.010x | 1.010x |
| TFG-CARMA3 | 25 | 1.006x | 1.003x | 1.081x |
| TFG-CARMA3 | 50 | 1.025x | 0.987x | 0.969x |
| TFG-CARMA3 | 100 | 1.040x | 0.989x | 0.961x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
