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

Source commit: `fee96e2b3be2d3d9d67b654f76ebea92b291df57`.

## CARMA3

CG baseline: 3-D displacement 0.2277 m, max roll/pitch 0.1177 deg, yaw 0.7469 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 0.999x | 1.020x | 1.004x | 1.001x | 1.010x | 0.0269 | 0.0269 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0269 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 0.998x | 1.002x | 1.006x | 1.009x | 1.003x | 0.0269 | n/a | n/a |
| mekf | x-forward | 10 | 0.997x | 1.004x | 1.006x | 1.011x | 1.002x | 0.0269 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 0.996x | 1.006x | 1.038x | 1.138x | 0.993x | 0.0269 | n/a | n/a |
| unmodeled | x-forward | 20 | 0.999x | 1.051x | 1.009x | 1.004x | 1.021x | 0.0538 | 0.0538 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0538 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 0.997x | 1.005x | 1.012x | 1.018x | 1.006x | 0.0538 | n/a | n/a |
| mekf | x-forward | 20 | 0.994x | 1.008x | 1.012x | 1.022x | 1.005x | 0.0538 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 0.993x | 1.012x | 1.103x | 1.300x | 0.986x | 0.0538 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.001x | 1.092x | 1.013x | 1.005x | 1.032x | 0.0807 | 0.0807 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0807 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 0.995x | 1.008x | 1.018x | 1.028x | 1.009x | 0.0807 | n/a | n/a |
| mekf | x-forward | 30 | 0.992x | 1.012x | 1.018x | 1.034x | 1.008x | 0.0807 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 0.991x | 1.020x | 1.176x | 1.482x | 0.977x | 0.0807 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.003x | 1.202x | 1.020x | 1.391x | 1.055x | 0.1346 | 0.1346 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1346 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.986x | 1.013x | 1.028x | 1.048x | 1.016x | 0.1346 | n/a | n/a |
| mekf | x-forward | 50 | 0.981x | 1.022x | 1.030x | 1.060x | 1.015x | 0.1346 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 0.981x | 1.041x | 1.346x | 1.896x | 0.984x | 0.1346 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.053x | 1.677x | 1.398x | 3.593x | 1.121x | 0.2691 | 0.2691 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2691 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 0.982x | 1.030x | 1.058x | 1.098x | 1.042x | 0.2691 | n/a | 135% |
| mekf | x-forward | 100 | 0.972x | 1.047x | 1.062x | 1.127x | 1.040x | 0.2691 | n/a | 154% |
| mekf_ou3 | x-forward | 100 | 0.983x | 1.145x | 2.045x | 3.494x | 2.785x | 0.2691 | n/a | 132% |
| unmodeled | x-forward | 200 | 1.252x | 3.271x | 3.247x | 8.106x | 1.300x | 0.5383 | 0.5383 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5383 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 0.998x | 1.093x | 1.026x | 1.071x | 1.082x | 0.5383 | n/a | 101% |
| mekf | x-forward | 200 | 0.972x | 1.133x | 1.143x | 1.297x | 1.126x | 0.5383 | n/a | 111% |
| mekf_ou3 | x-forward | 200 | 1.089x | 1.950x | 5.695x | 9.334x | 10.658x | 0.5383 | n/a | 65% |
| unmodeled | y-starboard | 10 | 1.002x | 1.005x | 1.121x | 1.157x | 1.043x | 0.0386 | 0.0386 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.002x | 1.004x | 1.003x | 1.003x | 1.002x | 0.0386 | n/a | n/a |
| mekf | y-starboard | 10 | 1.003x | 1.005x | 1.002x | 1.002x | 1.002x | 0.0386 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.002x | 1.005x | 0.985x | 1.010x | 0.979x | 0.0386 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.007x | 1.013x | 1.269x | 1.516x | 1.090x | 0.0772 | 0.0772 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.004x | 1.008x | 1.005x | 1.006x | 1.005x | 0.0772 | n/a | n/a |
| mekf | y-starboard | 20 | 1.007x | 1.012x | 1.003x | 1.003x | 1.004x | 0.0772 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.006x | 1.011x | 0.990x | 1.023x | 0.963x | 0.0772 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.012x | 1.023x | 1.435x | 2.006x | 1.141x | 0.1159 | 0.1159 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1159 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.006x | 1.010x | 1.007x | 1.007x | 1.008x | 0.1159 | n/a | n/a |
| mekf | y-starboard | 30 | 1.010x | 1.015x | 1.004x | 1.005x | 1.006x | 0.1159 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.009x | 1.014x | 0.995x | 1.040x | 0.951x | 0.1159 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.032x | 1.123x | 1.799x | 3.008x | 1.252x | 0.1931 | 0.1931 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1931 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.011x | 1.018x | 1.011x | 1.012x | 1.014x | 0.1931 | n/a | n/a |
| mekf | y-starboard | 50 | 1.018x | 1.027x | 1.007x | 1.011x | 1.013x | 0.1931 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.019x | 1.029x | 1.012x | 1.093x | 0.970x | 0.1931 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.115x | 1.737x | 2.785x | 5.508x | 1.562x | 0.3862 | 0.3862 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3862 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.023x | 1.042x | 1.020x | 1.027x | 1.038x | 0.3862 | n/a | 80% |
| mekf | y-starboard | 100 | 1.041x | 1.058x | 1.011x | 1.020x | 1.030x | 0.3862 | n/a | 65% |
| mekf_ou3 | y-starboard | 100 | 1.054x | 1.103x | 1.225x | 1.402x | 1.457x | 0.3862 | n/a | 53% |
| unmodeled | y-starboard | 200 | 1.418x | 3.609x | 4.774x | 10.253x | 2.233x | 0.7724 | 0.7724 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7724 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.088x | 1.127x | 1.006x | 1.099x | 1.040x | 0.7724 | n/a | 79% |
| mekf | y-starboard | 200 | 1.130x | 1.178x | 1.023x | 1.028x | 1.081x | 0.7724 | n/a | 69% |
| mekf_ou3 | y-starboard | 200 | 1.210x | 1.391x | 2.713x | 4.183x | 2.491x | 0.7724 | n/a | 50% |
| unmodeled | z-up | 10 | 1.001x | 1.043x | 1.002x | 1.021x | 1.000x | 0.0462 | 0.0462 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0462 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.003x | 1.000x | 1.002x | 1.003x | 0.0462 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.003x | 1.003x | 1.020x | 1.003x | 0.0462 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 1.003x | 1.004x | 1.017x | 1.002x | 0.0462 | n/a | n/a |
| unmodeled | z-up | 20 | 1.004x | 1.131x | 1.004x | 1.042x | 1.000x | 0.0924 | 0.0924 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0924 | 0.0000 | n/a |
| wrapper | z-up | 20 | 0.999x | 1.004x | 0.999x | 1.004x | 1.008x | 0.0924 | n/a | n/a |
| mekf | z-up | 20 | 0.999x | 1.004x | 1.005x | 1.040x | 1.008x | 0.0924 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.996x | 1.006x | 1.007x | 1.036x | 1.005x | 0.0924 | n/a | n/a |
| unmodeled | z-up | 30 | 1.009x | 1.253x | 1.006x | 1.063x | 0.999x | 0.1386 | 0.1386 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1386 | 0.0000 | n/a |
| wrapper | z-up | 30 | 0.998x | 1.006x | 0.998x | 1.004x | 1.015x | 0.1386 | n/a | n/a |
| mekf | z-up | 30 | 0.998x | 1.006x | 1.007x | 1.061x | 1.014x | 0.1386 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.994x | 1.011x | 1.011x | 1.055x | 1.009x | 0.1386 | n/a | n/a |
| unmodeled | z-up | 50 | 1.030x | 1.568x | 1.010x | 1.106x | 0.999x | 0.2311 | 0.2311 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2311 | 0.0000 | n/a |
| wrapper | z-up | 50 | 0.996x | 1.012x | 0.993x | 1.006x | 1.033x | 0.2311 | n/a | n/a |
| mekf | z-up | 50 | 0.996x | 1.012x | 1.009x | 1.103x | 1.032x | 0.2311 | n/a | n/a |
| mekf_ou3 | z-up | 50 | 0.991x | 1.023x | 1.017x | 1.096x | 1.019x | 0.2311 | n/a | n/a |
| unmodeled | z-up | 100 | 1.143x | 2.542x | 1.065x | 1.215x | 0.999x | 0.4621 | 0.4621 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.4621 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.012x | 1.041x | 0.982x | 1.009x | 1.112x | 0.4621 | n/a | 91% |
| mekf | z-up | 100 | 1.013x | 1.042x | 1.069x | 1.214x | 1.112x | 0.4621 | n/a | 91% |
| mekf_ou3 | z-up | 100 | 1.005x | 1.071x | 1.056x | 1.191x | 1.060x | 0.4621 | n/a | 97% |
| unmodeled | z-up | 200 | 1.508x | 4.700x | 1.332x | 2.013x | 1.003x | 0.9242 | 0.9242 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.9242 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.071x | 1.126x | 1.029x | 1.098x | 1.199x | 0.9242 | n/a | 86% |
| mekf | z-up | 200 | 1.065x | 1.122x | 1.353x | 2.029x | 1.475x | 0.9242 | n/a | 87% |
| mekf_ou3 | z-up | 200 | 1.060x | 1.227x | 1.326x | 2.029x | 1.230x | 0.9242 | n/a | 88% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2235 m, max roll/pitch 0.1157 deg, yaw 0.7324 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 0.999x | 1.020x | 1.007x | 1.000x | 1.011x | 0.0269 | 0.0269 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0269 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 0.998x | 1.003x | 1.003x | 1.004x | 1.001x | 0.0269 | n/a | n/a |
| mekf | x-forward | 10 | 0.997x | 1.004x | 1.003x | 1.006x | 1.001x | 0.0269 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 0.996x | 1.005x | 1.049x | 1.138x | 0.993x | 0.0269 | n/a | n/a |
| unmodeled | x-forward | 20 | 0.998x | 1.051x | 1.014x | 1.002x | 1.022x | 0.0538 | 0.0538 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0538 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 0.996x | 1.005x | 1.006x | 1.008x | 1.003x | 0.0538 | n/a | n/a |
| mekf | x-forward | 20 | 0.994x | 1.008x | 1.007x | 1.012x | 1.002x | 0.0538 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 0.993x | 1.012x | 1.115x | 1.302x | 0.986x | 0.0538 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.001x | 1.093x | 1.022x | 1.002x | 1.033x | 0.0807 | 0.0807 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0807 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 0.995x | 1.008x | 1.009x | 1.011x | 1.004x | 0.0807 | n/a | n/a |
| mekf | x-forward | 30 | 0.991x | 1.012x | 1.010x | 1.018x | 1.003x | 0.0807 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 0.991x | 1.020x | 1.189x | 1.484x | 0.978x | 0.0807 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.004x | 1.206x | 1.036x | 1.482x | 1.057x | 0.1346 | 0.1346 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1346 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.986x | 1.014x | 1.016x | 1.020x | 1.007x | 0.1346 | n/a | n/a |
| mekf | x-forward | 50 | 0.981x | 1.022x | 1.017x | 1.032x | 1.006x | 0.1346 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 0.982x | 1.041x | 1.361x | 1.896x | 0.992x | 0.1346 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.055x | 1.676x | 1.438x | 3.772x | 1.124x | 0.2691 | 0.2691 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2691 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 0.981x | 1.032x | 1.033x | 1.040x | 1.017x | 0.2691 | n/a | 135% |
| mekf | x-forward | 100 | 0.971x | 1.048x | 1.037x | 1.066x | 1.015x | 0.2691 | n/a | 152% |
| mekf_ou3 | x-forward | 100 | 0.985x | 1.130x | 1.909x | 3.122x | 1.698x | 0.2691 | n/a | 128% |
| unmodeled | x-forward | 200 | 1.261x | 3.270x | 3.352x | 8.468x | 1.296x | 0.5383 | 0.5383 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5383 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.001x | 1.095x | 1.007x | 1.059x | 1.019x | 0.5383 | n/a | 100% |
| mekf | x-forward | 200 | 0.972x | 1.132x | 1.082x | 1.149x | 1.039x | 0.5383 | n/a | 111% |
| mekf_ou3 | x-forward | 200 | 1.079x | 2.952x | 4.784x | 10.927x | 33.972x | 0.5383 | n/a | 70% |
| unmodeled | y-starboard | 10 | 1.002x | 1.006x | 1.131x | 1.155x | 1.047x | 0.0386 | 0.0386 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.002x | 1.005x | 1.002x | 1.002x | 1.001x | 0.0386 | n/a | n/a |
| mekf | y-starboard | 10 | 1.003x | 1.006x | 1.001x | 1.000x | 1.001x | 0.0386 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.003x | 1.005x | 0.995x | 1.009x | 0.985x | 0.0386 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.007x | 1.013x | 1.290x | 1.546x | 1.098x | 0.0772 | 0.0772 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.005x | 1.010x | 1.003x | 1.004x | 1.002x | 0.0772 | n/a | n/a |
| mekf | y-starboard | 20 | 1.007x | 1.012x | 1.001x | 1.001x | 1.001x | 0.0772 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.006x | 1.011x | 0.999x | 1.021x | 0.977x | 0.0772 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.013x | 1.025x | 1.468x | 2.060x | 1.152x | 0.1159 | 0.1159 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1159 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.007x | 1.012x | 1.004x | 1.005x | 1.004x | 0.1159 | n/a | n/a |
| mekf | y-starboard | 30 | 1.010x | 1.016x | 1.001x | 1.001x | 1.002x | 0.1159 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.009x | 1.015x | 1.005x | 1.039x | 0.978x | 0.1159 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.033x | 1.127x | 1.856x | 3.110x | 1.272x | 0.1931 | 0.1931 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1931 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.013x | 1.021x | 1.007x | 1.005x | 1.006x | 0.1931 | n/a | n/a |
| mekf | y-starboard | 50 | 1.019x | 1.029x | 1.002x | 1.001x | 1.005x | 0.1931 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.020x | 1.029x | 1.023x | 1.096x | 1.029x | 0.1931 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.120x | 1.737x | 2.900x | 5.734x | 1.605x | 0.3862 | 0.3862 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3862 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.029x | 1.044x | 1.014x | 1.012x | 1.016x | 0.3862 | n/a | 76% |
| mekf | y-starboard | 100 | 1.045x | 1.063x | 1.003x | 1.000x | 1.010x | 0.3862 | n/a | 63% |
| mekf_ou3 | y-starboard | 100 | 1.057x | 1.100x | 1.176x | 1.314x | 1.134x | 0.3862 | n/a | 52% |
| unmodeled | y-starboard | 200 | 1.429x | 3.610x | 4.997x | 10.723x | 2.320x | 0.7724 | 0.7724 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7724 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.102x | 1.140x | 1.028x | 1.108x | 1.000x | 0.7724 | n/a | 76% |
| mekf | y-starboard | 200 | 1.135x | 1.184x | 1.008x | 0.993x | 1.029x | 0.7724 | n/a | 69% |
| mekf_ou3 | y-starboard | 200 | 1.201x | 1.338x | 1.864x | 3.154x | 3.249x | 0.7724 | n/a | 53% |
| unmodeled | z-up | 10 | 1.002x | 1.043x | 1.002x | 1.020x | 1.000x | 0.0462 | 0.0462 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0462 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.001x | 1.004x | 0.999x | 1.002x | 1.002x | 0.0462 | n/a | n/a |
| mekf | z-up | 10 | 1.001x | 1.004x | 1.002x | 1.020x | 1.001x | 0.0462 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 1.000x | 1.004x | 1.003x | 1.019x | 1.001x | 0.0462 | n/a | n/a |
| unmodeled | z-up | 20 | 1.004x | 1.131x | 1.004x | 1.041x | 1.000x | 0.0924 | 0.0924 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0924 | 0.0000 | n/a |
| wrapper | z-up | 20 | 0.999x | 1.004x | 0.998x | 1.003x | 1.004x | 0.0924 | n/a | n/a |
| mekf | z-up | 20 | 0.999x | 1.004x | 1.003x | 1.042x | 1.003x | 0.0924 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.997x | 1.007x | 1.005x | 1.038x | 1.003x | 0.0924 | n/a | n/a |
| unmodeled | z-up | 30 | 1.011x | 1.253x | 1.007x | 1.060x | 1.000x | 0.1386 | 0.1386 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1386 | 0.0000 | n/a |
| wrapper | z-up | 30 | 0.998x | 1.006x | 0.996x | 1.002x | 1.006x | 0.1386 | n/a | n/a |
| mekf | z-up | 30 | 0.998x | 1.007x | 1.006x | 1.062x | 1.006x | 0.1386 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.995x | 1.011x | 1.009x | 1.059x | 1.005x | 0.1386 | n/a | n/a |
| unmodeled | z-up | 50 | 1.033x | 1.569x | 1.016x | 1.116x | 1.000x | 0.2311 | 0.2311 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2311 | 0.0000 | n/a |
| wrapper | z-up | 50 | 0.995x | 1.012x | 0.992x | 1.004x | 1.014x | 0.2311 | n/a | n/a |
| mekf | z-up | 50 | 0.998x | 1.012x | 1.013x | 1.115x | 1.015x | 0.2311 | n/a | n/a |
| mekf_ou3 | z-up | 50 | 0.993x | 1.024x | 1.019x | 1.120x | 1.011x | 0.2311 | n/a | n/a |
| unmodeled | z-up | 100 | 1.151x | 2.543x | 1.071x | 1.211x | 1.000x | 0.4621 | 0.4621 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.4621 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.013x | 1.043x | 0.990x | 1.002x | 1.045x | 0.4621 | n/a | 91% |
| mekf | z-up | 100 | 1.015x | 1.046x | 1.071x | 1.203x | 1.048x | 0.4621 | n/a | 90% |
| mekf_ou3 | z-up | 100 | 1.009x | 1.069x | 1.061x | 1.219x | 1.031x | 0.4621 | n/a | 94% |
| unmodeled | z-up | 200 | 1.520x | 4.711x | 1.339x | 2.081x | 1.004x | 0.9242 | 0.9242 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.9242 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.076x | 1.135x | 1.035x | 1.093x | 1.059x | 0.9242 | n/a | 85% |
| mekf | z-up | 200 | 1.065x | 1.124x | 1.344x | 2.082x | 1.196x | 0.9242 | n/a | 87% |
| mekf_ou3 | z-up | 200 | 1.066x | 1.222x | 1.328x | 2.088x | 1.114x | 0.9242 | n/a | 87% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.115x | 0.965x | 0.996x |
| CARMA3 | 2 | 1.035x | 0.974x | 0.998x |
| CARMA3 | 5 | 1.015x | 0.985x | 1.003x |
| CARMA3 | 10 | 1.013x | 0.986x | 1.032x |
| CARMA3 | 15 | 1.012x | 0.982x | 1.112x |
| CARMA3 | 25 | 1.023x | 1.035x | 1.197x |
| CARMA3 | 50 | 1.147x | 1.160x | 0.916x |
| CARMA3 | 100 | 1.140x | 1.208x | 0.951x |
| TFG-CARMA3 | 1 | 1.117x | 0.966x | 0.996x |
| TFG-CARMA3 | 2 | 1.035x | 0.969x | 0.996x |
| TFG-CARMA3 | 5 | 1.015x | 0.982x | 0.999x |
| TFG-CARMA3 | 10 | 1.013x | 0.987x | 1.012x |
| TFG-CARMA3 | 15 | 1.013x | 0.990x | 1.045x |
| TFG-CARMA3 | 25 | 1.027x | 1.043x | 1.060x |
| TFG-CARMA3 | 50 | 1.161x | 1.183x | 0.918x |
| TFG-CARMA3 | 100 | 1.154x | 1.228x | 0.954x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
