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

Source commit: `f8618df08cd9c2854dd1632be1a7a3728a562863`.

## CARMA3

CG baseline: 3-D displacement 0.2117 m, max roll/pitch 0.2539 deg, yaw 0.5805 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.001x | 1.103x | 1.005x | 1.013x | 0.999x | 0.0258 | 0.0258 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0258 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.002x | 0.999x | 1.002x | 1.001x | 0.0258 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.019x | 1.015x | 1.036x | 0.999x | 0.0258 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.022x | 1.075x | 1.142x | 0.989x | 0.0258 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.218x | 1.010x | 1.024x | 0.997x | 0.0516 | 0.0516 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0516 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.003x | 0.997x | 1.004x | 1.002x | 0.0516 | n/a | n/a |
| mekf | x-forward | 20 | 1.001x | 1.028x | 1.031x | 1.073x | 0.998x | 0.0516 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.034x | 1.157x | 1.296x | 0.978x | 0.0516 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.009x | 1.354x | 1.010x | 1.027x | 0.997x | 0.0773 | 0.0773 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0773 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.005x | 0.996x | 1.005x | 1.004x | 0.0773 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.042x | 1.042x | 1.098x | 0.999x | 0.0773 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.007x | 1.053x | 1.240x | 1.446x | 0.969x | 0.0773 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.026x | 1.648x | 1.012x | 1.046x | 0.995x | 0.1289 | 0.1289 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1289 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.007x | 0.994x | 1.003x | 1.009x | 0.1289 | n/a | n/a |
| mekf | x-forward | 50 | 1.005x | 1.056x | 1.074x | 1.184x | 1.001x | 0.1289 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.016x | 1.078x | 1.436x | 1.798x | 0.948x | 0.1289 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.093x | 2.499x | 1.006x | 1.096x | 1.006x | 0.2578 | 0.2578 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2578 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 1.000x | 1.013x | 0.991x | 0.997x | 1.032x | 0.2578 | n/a | 100% |
| mekf | x-forward | 100 | 1.013x | 1.151x | 1.169x | 1.380x | 1.031x | 0.2578 | n/a | 86% |
| mekf_ou3 | x-forward | 100 | 1.053x | 1.216x | 2.017x | 3.428x | 2.038x | 0.2578 | n/a | 43% |
| unmodeled | x-forward | 200 | 1.321x | 4.266x | 0.938x | 1.204x | 1.108x | 0.5155 | 0.5155 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5155 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.011x | 1.066x | 0.993x | 1.053x | 1.056x | 0.5155 | n/a | 97% |
| mekf | x-forward | 200 | 1.036x | 1.391x | 1.340x | 1.769x | 1.197x | 0.5155 | n/a | 89% |
| mekf_ou3 | x-forward | 200 | 1.178x | 1.613x | 4.075x | 10.443x | 7.685x | 0.5155 | n/a | 45% |
| unmodeled | y-starboard | 10 | 1.004x | 1.116x | 0.995x | 1.009x | 1.008x | 0.0298 | 0.0298 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0298 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.001x | 1.000x | 1.005x | 0.998x | 0.0298 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.006x | 1.002x | 1.009x | 0.994x | 0.0298 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.000x | 1.005x | 1.011x | 1.027x | 0.948x | 0.0298 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.011x | 1.256x | 0.991x | 1.018x | 1.017x | 0.0597 | 0.0597 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0597 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.001x | 1.001x | 1.001x | 1.010x | 0.995x | 0.0597 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.009x | 1.003x | 1.020x | 0.989x | 0.0597 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.001x | 1.010x | 1.023x | 1.054x | 0.923x | 0.0597 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.019x | 1.420x | 0.985x | 1.028x | 1.026x | 0.0895 | 0.0895 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0895 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.002x | 1.002x | 1.001x | 1.014x | 0.994x | 0.0895 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.006x | 1.004x | 1.026x | 0.984x | 0.0895 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.001x | 1.008x | 1.034x | 1.080x | 0.931x | 0.0895 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.043x | 1.793x | 0.972x | 1.046x | 1.051x | 0.1492 | 0.1492 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1492 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.003x | 1.003x | 1.003x | 1.024x | 0.993x | 0.1492 | n/a | n/a |
| mekf | y-starboard | 50 | 1.001x | 1.003x | 1.004x | 1.031x | 0.980x | 0.1492 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.007x | 1.008x | 1.058x | 1.131x | 1.062x | 0.1492 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.137x | 2.811x | 0.934x | 1.094x | 1.129x | 0.2984 | 0.2984 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2984 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.007x | 1.008x | 1.000x | 1.044x | 1.016x | 0.2984 | n/a | 95% |
| mekf | y-starboard | 100 | 1.007x | 1.018x | 0.996x | 1.035x | 0.994x | 0.2984 | n/a | 95% |
| mekf_ou3 | y-starboard | 100 | 1.036x | 1.044x | 1.276x | 2.249x | 2.525x | 0.2984 | n/a | 74% |
| unmodeled | y-starboard | 200 | 1.426x | 4.972x | 0.970x | 1.191x | 1.227x | 0.5968 | 0.5968 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5968 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.026x | 1.030x | 0.994x | 1.020x | 1.052x | 0.5968 | n/a | 94% |
| mekf | y-starboard | 200 | 1.034x | 1.041x | 0.935x | 1.027x | 1.135x | 0.5968 | n/a | 92% |
| mekf_ou3 | y-starboard | 200 | 1.177x | 1.182x | 4.268x | 6.984x | 8.538x | 0.5968 | n/a | 58% |
| unmodeled | z-up | 10 | 1.005x | 1.072x | 0.991x | 1.029x | 0.993x | 0.0386 | 0.0386 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.001x | 1.002x | 1.004x | 0.993x | 0.0386 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.002x | 0.994x | 1.032x | 0.986x | 0.0386 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 1.003x | 0.994x | 1.032x | 0.987x | 0.0386 | n/a | n/a |
| unmodeled | z-up | 20 | 1.015x | 1.200x | 0.982x | 1.061x | 0.987x | 0.0772 | 0.0772 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.002x | 1.003x | 1.008x | 0.987x | 0.0772 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.004x | 0.988x | 1.064x | 0.973x | 0.0772 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.998x | 1.007x | 0.988x | 1.066x | 0.974x | 0.0772 | n/a | n/a |
| unmodeled | z-up | 30 | 1.028x | 1.370x | 0.975x | 1.109x | 0.982x | 0.1158 | 0.1158 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1158 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.001x | 1.004x | 1.005x | 1.012x | 0.983x | 0.1158 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.005x | 0.984x | 1.114x | 0.963x | 0.1158 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.011x | 0.984x | 1.117x | 0.964x | 0.1158 | n/a | n/a |
| unmodeled | z-up | 50 | 1.064x | 1.784x | 0.959x | 1.205x | 0.977x | 0.1930 | 0.1930 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1930 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.001x | 1.007x | 1.009x | 1.021x | 0.979x | 0.1930 | n/a | 98% |
| mekf | z-up | 50 | 1.002x | 1.008x | 0.975x | 1.216x | 0.956x | 0.1930 | n/a | 97% |
| mekf_ou3 | z-up | 50 | 0.997x | 1.022x | 0.974x | 1.220x | 0.948x | 0.1930 | n/a | 105% |
| unmodeled | z-up | 100 | 1.208x | 2.997x | 0.926x | 1.452x | 0.994x | 0.3860 | 0.3860 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3860 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.019x | 1.017x | 1.038x | 1.036x | 0.3860 | n/a | 99% |
| mekf | z-up | 100 | 1.004x | 1.019x | 0.958x | 1.477x | 1.044x | 0.3860 | n/a | 98% |
| mekf_ou3 | z-up | 100 | 0.996x | 1.062x | 0.955x | 1.485x | 0.955x | 0.3860 | n/a | 102% |
| unmodeled | z-up | 200 | 1.640x | 5.605x | 0.892x | 1.953x | 1.166x | 0.7720 | 0.7720 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7720 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.014x | 1.070x | 1.030x | 1.081x | 1.181x | 0.7720 | n/a | 98% |
| mekf | z-up | 200 | 1.011x | 1.059x | 0.952x | 2.003x | 1.774x | 0.7720 | n/a | 98% |
| mekf_ou3 | z-up | 200 | 1.005x | 1.197x | 0.954x | 2.028x | 1.235x | 0.7720 | n/a | 99% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2118 m, max roll/pitch 0.2538 deg, yaw 0.5808 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.001x | 1.103x | 1.005x | 1.013x | 0.999x | 0.0258 | 0.0258 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0258 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.002x | 0.999x | 1.002x | 1.000x | 0.0258 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.019x | 1.015x | 1.036x | 0.999x | 0.0258 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.022x | 1.075x | 1.144x | 0.992x | 0.0258 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.218x | 1.010x | 1.024x | 0.998x | 0.0516 | 0.0516 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0516 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.003x | 0.998x | 1.005x | 1.001x | 0.0516 | n/a | n/a |
| mekf | x-forward | 20 | 1.001x | 1.028x | 1.031x | 1.073x | 0.997x | 0.0516 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.034x | 1.158x | 1.299x | 0.988x | 0.0516 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.009x | 1.354x | 1.010x | 1.026x | 0.998x | 0.0773 | 0.0773 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0773 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 0.999x | 1.005x | 0.997x | 1.007x | 1.001x | 0.0773 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.042x | 1.042x | 1.098x | 0.997x | 0.0773 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.007x | 1.053x | 1.242x | 1.451x | 0.991x | 0.0773 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.026x | 1.648x | 1.012x | 1.045x | 0.996x | 0.1289 | 0.1289 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1289 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.007x | 0.995x | 1.010x | 1.003x | 0.1289 | n/a | n/a |
| mekf | x-forward | 50 | 1.004x | 1.056x | 1.074x | 1.184x | 0.994x | 0.1289 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.017x | 1.078x | 1.436x | 1.805x | 1.007x | 0.1289 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.093x | 2.499x | 1.005x | 1.094x | 1.007x | 0.2578 | 0.2578 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2578 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 1.000x | 1.013x | 0.991x | 1.017x | 1.008x | 0.2578 | n/a | 100% |
| mekf | x-forward | 100 | 1.012x | 1.151x | 1.168x | 1.382x | 1.004x | 0.2578 | n/a | 87% |
| mekf_ou3 | x-forward | 100 | 1.048x | 1.216x | 2.016x | 2.810x | 1.280x | 0.2578 | n/a | 49% |
| unmodeled | x-forward | 200 | 1.321x | 4.267x | 0.938x | 1.203x | 1.111x | 0.5155 | 0.5155 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5155 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.011x | 1.066x | 0.991x | 1.108x | 1.009x | 0.5155 | n/a | 97% |
| mekf | x-forward | 200 | 1.035x | 1.391x | 1.339x | 1.777x | 1.082x | 0.5155 | n/a | 89% |
| mekf_ou3 | x-forward | 200 | 1.181x | 1.544x | 2.511x | 3.486x | 11.876x | 0.5155 | n/a | 44% |
| unmodeled | y-starboard | 10 | 1.004x | 1.116x | 0.995x | 1.009x | 1.008x | 0.0298 | 0.0298 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0298 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.001x | 1.000x | 1.003x | 0.999x | 0.0298 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.006x | 1.001x | 1.009x | 0.996x | 0.0298 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.000x | 1.005x | 1.011x | 1.028x | 0.947x | 0.0298 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.011x | 1.256x | 0.991x | 1.018x | 1.017x | 0.0597 | 0.0597 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0597 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.001x | 1.001x | 1.000x | 1.006x | 0.998x | 0.0597 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.009x | 1.002x | 1.020x | 0.992x | 0.0597 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.001x | 1.010x | 1.022x | 1.055x | 0.921x | 0.0597 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.019x | 1.420x | 0.985x | 1.028x | 1.027x | 0.0895 | 0.0895 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0895 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.001x | 1.002x | 0.999x | 1.009x | 0.997x | 0.0895 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.006x | 1.002x | 1.027x | 0.988x | 0.0895 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.001x | 1.008x | 1.033x | 1.080x | 0.927x | 0.0895 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.043x | 1.793x | 0.972x | 1.047x | 1.053x | 0.1492 | 0.1492 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1492 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.003x | 1.003x | 1.000x | 1.014x | 0.996x | 0.1492 | n/a | n/a |
| mekf | y-starboard | 50 | 1.001x | 1.003x | 1.002x | 1.032x | 0.982x | 0.1492 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.006x | 1.007x | 1.057x | 1.133x | 1.038x | 0.1492 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.137x | 2.811x | 0.934x | 1.105x | 1.133x | 0.2984 | 0.2984 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2984 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.007x | 1.008x | 0.995x | 1.026x | 0.998x | 0.2984 | n/a | 95% |
| mekf | y-starboard | 100 | 1.007x | 1.018x | 0.992x | 1.039x | 0.973x | 0.2984 | n/a | 95% |
| mekf_ou3 | y-starboard | 100 | 1.031x | 1.039x | 1.238x | 2.206x | 1.895x | 0.2984 | n/a | 77% |
| unmodeled | y-starboard | 200 | 1.426x | 4.973x | 0.976x | 1.206x | 1.230x | 0.5968 | 0.5968 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5968 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.025x | 1.029x | 0.992x | 1.029x | 1.003x | 0.5968 | n/a | 94% |
| mekf | y-starboard | 200 | 1.032x | 1.039x | 0.933x | 0.998x | 1.001x | 0.5968 | n/a | 93% |
| mekf_ou3 | y-starboard | 200 | 1.134x | 1.180x | 2.587x | 5.125x | 8.457x | 0.5968 | n/a | 69% |
| unmodeled | z-up | 10 | 1.005x | 1.072x | 0.992x | 1.029x | 0.993x | 0.0386 | 0.0386 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.001x | 1.001x | 1.002x | 0.996x | 0.0386 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.002x | 0.993x | 1.030x | 0.989x | 0.0386 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 1.003x | 0.993x | 1.031x | 0.990x | 0.0386 | n/a | n/a |
| unmodeled | z-up | 20 | 1.015x | 1.200x | 0.984x | 1.061x | 0.986x | 0.0772 | 0.0772 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.002x | 1.002x | 1.004x | 0.993x | 0.0772 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.004x | 0.987x | 1.062x | 0.978x | 0.0772 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.998x | 1.007x | 0.987x | 1.065x | 0.981x | 0.0772 | n/a | n/a |
| unmodeled | z-up | 30 | 1.028x | 1.370x | 0.978x | 1.109x | 0.980x | 0.1158 | 0.1158 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1158 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.001x | 1.004x | 1.003x | 1.006x | 0.990x | 0.1158 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.005x | 0.982x | 1.111x | 0.968x | 0.1158 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.011x | 0.982x | 1.115x | 0.973x | 0.1158 | n/a | n/a |
| unmodeled | z-up | 50 | 1.064x | 1.784x | 0.963x | 1.205x | 0.973x | 0.1930 | 0.1930 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1930 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.001x | 1.007x | 1.005x | 1.011x | 0.984x | 0.1930 | n/a | 99% |
| mekf | z-up | 50 | 1.002x | 1.008x | 0.971x | 1.209x | 0.954x | 0.1930 | n/a | 97% |
| mekf_ou3 | z-up | 50 | 0.997x | 1.022x | 0.972x | 1.215x | 0.959x | 0.1930 | n/a | 105% |
| unmodeled | z-up | 100 | 1.208x | 2.998x | 0.934x | 1.451x | 0.975x | 0.3860 | 0.3860 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3860 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.019x | 1.010x | 1.021x | 0.990x | 0.3860 | n/a | 99% |
| mekf | z-up | 100 | 1.004x | 1.020x | 0.950x | 1.459x | 0.964x | 0.3860 | n/a | 98% |
| mekf_ou3 | z-up | 100 | 0.997x | 1.062x | 0.951x | 1.474x | 0.947x | 0.3860 | n/a | 102% |
| unmodeled | z-up | 200 | 1.640x | 5.608x | 0.911x | 1.952x | 1.079x | 0.7720 | 0.7720 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7720 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.014x | 1.070x | 1.022x | 1.065x | 1.049x | 0.7720 | n/a | 98% |
| mekf | z-up | 200 | 1.011x | 1.059x | 0.941x | 1.970x | 1.294x | 0.7720 | n/a | 98% |
| mekf_ou3 | z-up | 200 | 1.006x | 1.197x | 0.947x | 2.003x | 1.059x | 0.7720 | n/a | 99% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.123x | 1.032x | 1.030x |
| CARMA3 | 2 | 1.029x | 1.025x | 1.018x |
| CARMA3 | 5 | 1.006x | 1.014x | 0.999x |
| CARMA3 | 10 | 1.003x | 1.013x | 0.985x |
| CARMA3 | 15 | 1.002x | 1.017x | 1.036x |
| CARMA3 | 25 | 1.013x | 1.025x | 1.154x |
| CARMA3 | 50 | 1.092x | 1.032x | 1.007x |
| CARMA3 | 100 | 1.090x | 1.028x | 1.004x |
| TFG-CARMA3 | 1 | 1.123x | 1.035x | 1.024x |
| TFG-CARMA3 | 2 | 1.029x | 1.026x | 1.017x |
| TFG-CARMA3 | 5 | 1.006x | 1.013x | 1.003x |
| TFG-CARMA3 | 10 | 1.003x | 1.009x | 0.988x |
| TFG-CARMA3 | 15 | 1.002x | 1.010x | 0.990x |
| TFG-CARMA3 | 25 | 1.013x | 1.017x | 1.029x |
| TFG-CARMA3 | 50 | 1.092x | 1.031x | 0.998x |
| TFG-CARMA3 | 100 | 1.089x | 1.027x | 0.994x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
