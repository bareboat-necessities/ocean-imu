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

Source commit: `5b66d7aea7aaeaf8061e2208b904912cb2b65462`.

## CARMA3

CG baseline: 3-D displacement 0.2316 m, max roll/pitch 0.1858 deg, yaw 0.7131 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 0.999x | 1.063x | 0.990x | 0.998x | 1.000x | 0.0161 | 0.0161 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0161 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 0.998x | 1.002x | 1.008x | 0.0161 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.002x | 1.002x | 1.012x | 1.007x | 0.0161 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.005x | 1.103x | 1.152x | 1.000x | 0.0161 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.000x | 1.156x | 0.979x | 0.996x | 0.999x | 0.0322 | 0.0322 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0322 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.002x | 0.996x | 1.004x | 1.013x | 0.0322 | n/a | n/a |
| mekf | x-forward | 20 | 1.000x | 1.007x | 1.005x | 1.025x | 1.015x | 0.0322 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.013x | 1.215x | 1.318x | 1.010x | 0.0322 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.001x | 1.271x | 0.982x | 0.994x | 0.998x | 0.0483 | 0.0483 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0483 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.003x | 0.993x | 1.006x | 1.024x | 0.0483 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.014x | 1.007x | 1.039x | 1.024x | 0.0483 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.007x | 1.024x | 1.335x | 1.495x | 1.012x | 0.0483 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.008x | 1.541x | 0.996x | 1.031x | 0.997x | 0.0805 | 0.0805 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0805 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 1.000x | 1.006x | 0.990x | 1.009x | 1.044x | 0.0805 | n/a | n/a |
| mekf | x-forward | 50 | 1.001x | 1.032x | 1.014x | 1.067x | 1.049x | 0.0805 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.014x | 1.052x | 1.598x | 1.873x | 1.021x | 0.0805 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.047x | 2.289x | 1.044x | 1.179x | 0.998x | 0.1609 | 0.1609 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1609 | 0.0000 | n/a |
| wrapper | x-forward | 100 | 1.001x | 1.016x | 0.982x | 1.029x | 1.140x | 0.1609 | n/a | n/a |
| mekf | x-forward | 100 | 1.004x | 1.064x | 1.021x | 1.109x | 1.138x | 0.1609 | n/a | n/a |
| mekf_ou3 | x-forward | 100 | 1.048x | 1.138x | 2.825x | 4.182x | 2.552x | 0.1609 | n/a | n/a |
| unmodeled | x-forward | 200 | 1.200x | 3.806x | 1.077x | 1.265x | 0.999x | 0.3218 | 0.3218 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3218 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.008x | 1.061x | 0.981x | 1.054x | 1.206x | 0.3218 | n/a | 96% |
| mekf | x-forward | 200 | 1.018x | 1.060x | 1.037x | 1.200x | 1.435x | 0.3218 | n/a | 91% |
| mekf_ou3 | x-forward | 200 | 1.319x | 1.650x | 6.156x | 7.699x | 10.305x | 0.3218 | n/a | -60% |
| unmodeled | y-starboard | 10 | 1.004x | 1.011x | 0.997x | 1.009x | 1.003x | 0.0087 | 0.0087 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0087 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.004x | 1.006x | 0.991x | 0.0087 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.000x | 1.004x | 1.009x | 0.990x | 0.0087 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 0.999x | 1.017x | 1.019x | 0.967x | 0.0087 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.009x | 1.030x | 0.995x | 1.017x | 1.003x | 0.0173 | 0.0173 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0173 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.007x | 1.012x | 0.980x | 0.0173 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.000x | 1.007x | 1.014x | 0.981x | 0.0173 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.997x | 0.999x | 1.034x | 1.039x | 0.962x | 0.0173 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.016x | 1.068x | 1.002x | 1.025x | 1.003x | 0.0260 | 0.0260 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0260 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.009x | 1.017x | 0.975x | 0.0260 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.000x | 1.009x | 1.019x | 0.973x | 0.0260 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.997x | 0.999x | 1.052x | 1.060x | 0.984x | 0.0260 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.033x | 1.207x | 1.024x | 1.042x | 1.006x | 0.0434 | 0.0434 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0434 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 0.999x | 1.000x | 1.013x | 1.025x | 0.970x | 0.0434 | n/a | n/a |
| mekf | y-starboard | 50 | 1.000x | 1.000x | 1.013x | 1.029x | 0.969x | 0.0434 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.997x | 1.001x | 1.100x | 1.122x | 1.102x | 0.0434 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.093x | 1.747x | 1.072x | 1.096x | 1.011x | 0.0867 | 0.0867 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0867 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.000x | 1.002x | 1.011x | 1.030x | 0.989x | 0.0867 | n/a | 100% |
| mekf | y-starboard | 100 | 1.001x | 1.002x | 1.009x | 1.040x | 0.983x | 0.0867 | n/a | 99% |
| mekf_ou3 | y-starboard | 100 | 1.011x | 1.015x | 1.900x | 3.023x | 2.139x | 0.0867 | n/a | 88% |
| unmodeled | y-starboard | 200 | 1.279x | 3.105x | 1.140x | 1.191x | 1.025x | 0.1734 | 0.1734 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1734 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.005x | 1.031x | 0.988x | 1.002x | 1.072x | 0.1734 | n/a | 98% |
| mekf | y-starboard | 200 | 1.007x | 1.063x | 0.967x | 1.011x | 1.281x | 0.1734 | n/a | 98% |
| mekf_ou3 | y-starboard | 200 | 1.265x | 1.339x | 7.178x | 8.316x | 7.868x | 0.1734 | n/a | 5% |
| unmodeled | z-up | 10 | 1.011x | 1.041x | 1.001x | 1.004x | 1.000x | 0.0181 | 0.0181 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0181 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.004x | 1.006x | 0.989x | 0.0181 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.000x | 1.004x | 1.009x | 0.987x | 0.0181 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 0.999x | 1.004x | 1.007x | 0.991x | 0.0181 | n/a | n/a |
| unmodeled | z-up | 20 | 1.024x | 1.120x | 1.001x | 1.006x | 0.997x | 0.0362 | 0.0362 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0362 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.001x | 1.007x | 1.011x | 0.976x | 0.0362 | n/a | n/a |
| mekf | z-up | 20 | 1.000x | 1.001x | 1.008x | 1.017x | 0.974x | 0.0362 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.997x | 0.997x | 1.007x | 1.013x | 0.981x | 0.0362 | n/a | n/a |
| unmodeled | z-up | 30 | 1.039x | 1.230x | 1.001x | 1.009x | 0.997x | 0.0543 | 0.0543 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0543 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.000x | 1.001x | 1.011x | 1.018x | 0.969x | 0.0543 | n/a | n/a |
| mekf | z-up | 30 | 1.000x | 1.001x | 1.013x | 1.026x | 0.966x | 0.0543 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.996x | 0.996x | 1.011x | 1.020x | 0.973x | 0.0543 | n/a | n/a |
| unmodeled | z-up | 50 | 1.075x | 1.516x | 1.003x | 1.017x | 0.996x | 0.0904 | 0.0904 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0904 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.001x | 1.002x | 1.018x | 1.032x | 0.976x | 0.0904 | n/a | 99% |
| mekf | z-up | 50 | 1.001x | 1.002x | 1.021x | 1.047x | 0.973x | 0.0904 | n/a | 99% |
| mekf_ou3 | z-up | 50 | 0.993x | 0.995x | 1.019x | 1.036x | 0.962x | 0.0904 | n/a | 109% |
| unmodeled | z-up | 100 | 1.197x | 2.410x | 1.006x | 1.034x | 0.989x | 0.1808 | 0.1808 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1808 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.005x | 1.025x | 1.048x | 1.197x | 0.1808 | n/a | 99% |
| mekf | z-up | 100 | 1.002x | 1.006x | 1.040x | 1.084x | 1.130x | 0.1808 | n/a | 99% |
| mekf_ou3 | z-up | 100 | 0.987x | 0.998x | 1.039x | 1.073x | 0.970x | 0.1808 | n/a | 107% |
| unmodeled | z-up | 200 | 1.532x | 4.411x | 1.014x | 1.077x | 0.983x | 0.3617 | 0.3617 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3617 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.008x | 1.033x | 1.059x | 1.407x | 1.670x | 0.3617 | n/a | 98% |
| mekf | z-up | 200 | 1.014x | 1.027x | 1.360x | 2.114x | 2.997x | 0.3617 | n/a | 97% |
| mekf_ou3 | z-up | 200 | 0.979x | 1.026x | 1.063x | 1.129x | 1.410x | 0.3617 | n/a | 104% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2317 m, max roll/pitch 0.1862 deg, yaw 0.7095 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 0.999x | 1.063x | 0.990x | 0.998x | 0.998x | 0.0161 | 0.0161 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0161 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 0.999x | 1.001x | 1.003x | 0.0161 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.002x | 1.003x | 1.012x | 1.003x | 0.0161 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.005x | 1.105x | 1.154x | 0.998x | 0.0161 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.000x | 1.156x | 0.979x | 0.996x | 0.998x | 0.0322 | 0.0322 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0322 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.002x | 0.997x | 1.002x | 1.005x | 0.0322 | n/a | n/a |
| mekf | x-forward | 20 | 1.000x | 1.007x | 1.006x | 1.025x | 1.007x | 0.0322 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.013x | 1.220x | 1.324x | 1.006x | 0.0322 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.001x | 1.271x | 0.977x | 0.994x | 0.996x | 0.0483 | 0.0483 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0483 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.003x | 0.996x | 1.004x | 1.009x | 0.0483 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.014x | 1.010x | 1.039x | 1.011x | 0.0483 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.007x | 1.024x | 1.342x | 1.505x | 1.007x | 0.0483 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.008x | 1.541x | 0.992x | 1.033x | 0.995x | 0.0805 | 0.0805 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0805 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 1.000x | 1.006x | 0.993x | 1.007x | 1.017x | 0.0805 | n/a | n/a |
| mekf | x-forward | 50 | 1.001x | 1.032x | 1.017x | 1.067x | 1.021x | 0.0805 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.014x | 1.052x | 1.606x | 1.888x | 1.024x | 0.0805 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.047x | 2.289x | 1.040x | 1.181x | 0.999x | 0.1609 | 0.1609 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1609 | 0.0000 | n/a |
| wrapper | x-forward | 100 | 1.001x | 1.016x | 0.988x | 1.022x | 1.051x | 0.1609 | n/a | n/a |
| mekf | x-forward | 100 | 1.004x | 1.064x | 1.026x | 1.111x | 1.057x | 0.1609 | n/a | n/a |
| mekf_ou3 | x-forward | 100 | 1.041x | 1.126x | 2.537x | 3.167x | 1.434x | 0.1609 | n/a | n/a |
| unmodeled | x-forward | 200 | 1.200x | 3.806x | 1.074x | 1.267x | 0.997x | 0.3218 | 0.3218 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3218 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.007x | 1.062x | 0.983x | 1.048x | 1.071x | 0.3218 | n/a | 97% |
| mekf | x-forward | 200 | 1.014x | 1.060x | 0.999x | 1.086x | 1.168x | 0.3218 | n/a | 93% |
| mekf_ou3 | x-forward | 200 | 1.145x | 1.275x | 3.512x | 9.104x | 5.807x | 0.3218 | n/a | 27% |
| unmodeled | y-starboard | 10 | 1.004x | 1.012x | 0.997x | 1.009x | 1.002x | 0.0087 | 0.0087 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0087 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.002x | 1.003x | 0.994x | 0.0087 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.000x | 1.002x | 1.006x | 0.995x | 0.0087 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 0.999x | 1.016x | 1.019x | 0.974x | 0.0087 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.009x | 1.030x | 0.995x | 1.017x | 1.002x | 0.0173 | 0.0173 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0173 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.003x | 1.005x | 0.989x | 0.0173 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.000x | 1.003x | 1.007x | 0.989x | 0.0173 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.997x | 0.999x | 1.031x | 1.039x | 0.972x | 0.0173 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.016x | 1.068x | 0.998x | 1.025x | 1.002x | 0.0260 | 0.0260 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0260 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.004x | 1.008x | 0.986x | 0.0260 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.000x | 1.005x | 1.010x | 0.984x | 0.0260 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.997x | 0.999x | 1.047x | 1.060x | 0.996x | 0.0260 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.032x | 1.207x | 1.020x | 1.043x | 1.005x | 0.0434 | 0.0434 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0434 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.000x | 1.000x | 1.005x | 1.010x | 0.980x | 0.0434 | n/a | n/a |
| mekf | y-starboard | 50 | 1.000x | 1.000x | 1.006x | 1.015x | 0.979x | 0.0434 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.997x | 1.000x | 1.082x | 1.107x | 1.115x | 0.0434 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.093x | 1.746x | 1.069x | 1.100x | 1.013x | 0.0867 | 0.0867 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0867 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 0.999x | 1.002x | 1.003x | 1.010x | 0.978x | 0.0867 | n/a | 101% |
| mekf | y-starboard | 100 | 1.001x | 1.002x | 1.002x | 1.020x | 0.976x | 0.0867 | n/a | 99% |
| mekf_ou3 | y-starboard | 100 | 1.010x | 1.032x | 1.624x | 2.447x | 1.802x | 0.0867 | n/a | 89% |
| unmodeled | y-starboard | 200 | 1.279x | 3.105x | 1.137x | 1.195x | 1.025x | 0.1734 | 0.1734 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1734 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.004x | 1.031x | 0.987x | 1.002x | 0.993x | 0.1734 | n/a | 99% |
| mekf | y-starboard | 200 | 1.004x | 1.063x | 0.971x | 0.997x | 1.023x | 0.1734 | n/a | 99% |
| mekf_ou3 | y-starboard | 200 | 1.142x | 3.391x | 4.210x | 9.969x | 38.816x | 0.1734 | n/a | 49% |
| unmodeled | z-up | 10 | 1.011x | 1.041x | 1.001x | 1.004x | 0.999x | 0.0181 | 0.0181 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0181 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.003x | 0.994x | 0.0181 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.006x | 0.992x | 0.0181 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.998x | 0.999x | 1.003x | 1.006x | 0.994x | 0.0181 | n/a | n/a |
| unmodeled | z-up | 20 | 1.023x | 1.120x | 1.001x | 1.006x | 0.997x | 0.0362 | 0.0362 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0362 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.000x | 1.003x | 1.005x | 0.985x | 0.0362 | n/a | n/a |
| mekf | z-up | 20 | 1.000x | 1.001x | 1.004x | 1.011x | 0.983x | 0.0362 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.997x | 0.997x | 1.005x | 1.009x | 0.987x | 0.0362 | n/a | n/a |
| unmodeled | z-up | 30 | 1.038x | 1.230x | 1.002x | 1.009x | 0.995x | 0.0543 | 0.0543 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0543 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.000x | 1.001x | 1.004x | 1.007x | 0.976x | 0.0543 | n/a | n/a |
| mekf | z-up | 30 | 1.000x | 1.001x | 1.007x | 1.016x | 0.973x | 0.0543 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.996x | 0.996x | 1.007x | 1.014x | 0.980x | 0.0543 | n/a | n/a |
| unmodeled | z-up | 50 | 1.075x | 1.516x | 1.003x | 1.017x | 0.994x | 0.0904 | 0.0904 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0904 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.000x | 1.001x | 1.007x | 1.012x | 0.965x | 0.0904 | n/a | 99% |
| mekf | z-up | 50 | 1.000x | 1.002x | 1.011x | 1.029x | 0.960x | 0.0904 | n/a | 99% |
| mekf_ou3 | z-up | 50 | 0.993x | 0.995x | 1.013x | 1.026x | 0.968x | 0.0904 | n/a | 109% |
| unmodeled | z-up | 100 | 1.197x | 2.410x | 1.007x | 1.034x | 0.987x | 0.1808 | 0.1808 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1808 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.001x | 1.004x | 1.013x | 1.024x | 1.006x | 0.1808 | n/a | 99% |
| mekf | z-up | 100 | 1.002x | 1.005x | 1.022x | 1.058x | 0.988x | 0.1808 | n/a | 99% |
| mekf_ou3 | z-up | 100 | 0.987x | 0.998x | 1.026x | 1.051x | 0.942x | 0.1808 | n/a | 107% |
| unmodeled | z-up | 200 | 1.532x | 4.412x | 1.017x | 1.077x | 0.977x | 0.3617 | 0.3617 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3617 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.007x | 1.033x | 1.008x | 1.057x | 1.175x | 0.3617 | n/a | 99% |
| mekf | z-up | 200 | 1.007x | 1.023x | 1.032x | 1.134x | 1.714x | 0.3617 | n/a | 99% |
| mekf_ou3 | z-up | 200 | 0.978x | 1.026x | 1.053x | 1.109x | 1.070x | 0.3617 | n/a | 104% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.026x | 1.008x | 0.999x |
| CARMA3 | 2 | 1.007x | 1.004x | 0.997x |
| CARMA3 | 5 | 1.002x | 1.006x | 0.986x |
| CARMA3 | 10 | 1.001x | 1.018x | 0.985x |
| CARMA3 | 15 | 1.002x | 1.025x | 1.197x |
| CARMA3 | 25 | 1.009x | 1.023x | 1.768x |
| CARMA3 | 50 | 1.024x | 1.008x | 0.958x |
| CARMA3 | 100 | 1.024x | 1.017x | 0.945x |
| TFG-CARMA3 | 1 | 1.026x | 1.008x | 0.998x |
| TFG-CARMA3 | 2 | 1.007x | 1.004x | 0.997x |
| TFG-CARMA3 | 5 | 1.002x | 1.003x | 0.991x |
| TFG-CARMA3 | 10 | 1.001x | 1.008x | 0.968x |
| TFG-CARMA3 | 15 | 1.001x | 1.013x | 1.006x |
| TFG-CARMA3 | 25 | 1.007x | 1.008x | 1.124x |
| TFG-CARMA3 | 50 | 1.023x | 1.013x | 0.958x |
| TFG-CARMA3 | 100 | 1.024x | 1.017x | 0.949x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
