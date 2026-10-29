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

Source commit: `332964a53d21deae85257515bf54299e604aff06`.

## CARMA3

CG baseline: 3-D displacement 0.2124 m, max roll/pitch 0.2554 deg, yaw 0.5817 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.001x | 1.103x | 1.005x | 1.013x | 0.999x | 0.0258 | 0.0258 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0258 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.002x | 0.998x | 1.002x | 1.001x | 0.0258 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.018x | 1.013x | 1.034x | 0.999x | 0.0258 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.001x | 1.022x | 1.073x | 1.141x | 0.988x | 0.0258 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.218x | 1.010x | 1.023x | 0.997x | 0.0516 | 0.0516 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0516 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 0.999x | 1.003x | 0.995x | 1.004x | 1.003x | 0.0516 | n/a | n/a |
| mekf | x-forward | 20 | 1.000x | 1.027x | 1.028x | 1.068x | 0.999x | 0.0516 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.003x | 1.034x | 1.153x | 1.294x | 0.977x | 0.0516 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.009x | 1.354x | 1.010x | 1.025x | 0.997x | 0.0773 | 0.0773 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0773 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 0.999x | 1.005x | 0.992x | 1.005x | 1.004x | 0.0773 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.042x | 1.037x | 1.090x | 0.999x | 0.0773 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.006x | 1.053x | 1.235x | 1.442x | 0.967x | 0.0773 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.025x | 1.648x | 1.012x | 1.043x | 0.996x | 0.1289 | 0.1289 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1289 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.007x | 0.988x | 1.003x | 1.010x | 0.1289 | n/a | n/a |
| mekf | x-forward | 50 | 1.004x | 1.055x | 1.066x | 1.170x | 1.002x | 0.1289 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.015x | 1.078x | 1.426x | 1.791x | 0.945x | 0.1289 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.092x | 2.500x | 1.006x | 1.089x | 1.008x | 0.2578 | 0.2578 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2578 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 1.000x | 1.013x | 0.978x | 0.993x | 1.033x | 0.2578 | n/a | 100% |
| mekf | x-forward | 100 | 1.012x | 1.150x | 1.151x | 1.352x | 1.030x | 0.2578 | n/a | 87% |
| mekf_ou3 | x-forward | 100 | 1.049x | 1.215x | 1.996x | 3.296x | 2.015x | 0.2578 | n/a | 46% |
| unmodeled | x-forward | 200 | 1.318x | 4.269x | 0.939x | 1.188x | 1.111x | 0.5155 | 0.5155 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5155 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.017x | 1.068x | 0.975x | 1.049x | 1.058x | 0.5155 | n/a | 95% |
| mekf | x-forward | 200 | 1.035x | 1.388x | 1.305x | 1.709x | 1.189x | 0.5155 | n/a | 89% |
| mekf_ou3 | x-forward | 200 | 1.169x | 1.613x | 4.031x | 10.263x | 7.700x | 0.5155 | n/a | 47% |
| unmodeled | y-starboard | 10 | 1.004x | 1.116x | 0.995x | 1.009x | 1.008x | 0.0298 | 0.0298 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0298 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.001x | 1.000x | 1.003x | 0.998x | 0.0298 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.006x | 1.001x | 1.009x | 0.995x | 0.0298 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.000x | 1.005x | 1.011x | 1.027x | 0.946x | 0.0298 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.011x | 1.256x | 0.991x | 1.017x | 1.017x | 0.0597 | 0.0597 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0597 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.001x | 1.001x | 1.001x | 1.005x | 0.997x | 0.0597 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.009x | 1.003x | 1.019x | 0.991x | 0.0597 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.001x | 1.010x | 1.022x | 1.054x | 0.920x | 0.0597 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.019x | 1.421x | 0.985x | 1.026x | 1.027x | 0.0895 | 0.0895 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0895 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.001x | 1.002x | 1.000x | 1.007x | 0.996x | 0.0895 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.005x | 1.003x | 1.025x | 0.986x | 0.0895 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.002x | 1.008x | 1.033x | 1.079x | 0.925x | 0.0895 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.043x | 1.794x | 0.972x | 1.044x | 1.052x | 0.1492 | 0.1492 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1492 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.003x | 1.003x | 1.002x | 1.009x | 0.996x | 0.1492 | n/a | n/a |
| mekf | y-starboard | 50 | 1.001x | 1.003x | 1.003x | 1.029x | 0.982x | 0.1492 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.008x | 1.009x | 1.057x | 1.130x | 1.050x | 0.1492 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.137x | 2.813x | 0.935x | 1.099x | 1.131x | 0.2984 | 0.2984 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2984 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.007x | 1.008x | 0.998x | 1.007x | 1.020x | 0.2984 | n/a | 95% |
| mekf | y-starboard | 100 | 1.007x | 1.018x | 0.994x | 1.030x | 0.995x | 0.2984 | n/a | 95% |
| mekf_ou3 | y-starboard | 100 | 1.037x | 1.046x | 1.254x | 2.226x | 2.494x | 0.2984 | n/a | 73% |
| unmodeled | y-starboard | 200 | 1.425x | 4.976x | 0.983x | 1.193x | 1.227x | 0.5968 | 0.5968 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5968 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.033x | 1.038x | 0.992x | 1.017x | 1.059x | 0.5968 | n/a | 92% |
| mekf | y-starboard | 200 | 1.033x | 1.040x | 0.932x | 0.989x | 1.134x | 0.5968 | n/a | 92% |
| mekf_ou3 | y-starboard | 200 | 1.179x | 1.188x | 4.203x | 6.964x | 8.489x | 0.5968 | n/a | 58% |
| unmodeled | z-up | 10 | 1.005x | 1.072x | 0.991x | 1.009x | 0.993x | 0.0386 | 0.0386 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.001x | 1.002x | 1.004x | 0.993x | 0.0386 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.002x | 0.994x | 1.002x | 0.985x | 0.0386 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 1.003x | 0.993x | 1.006x | 0.987x | 0.0386 | n/a | n/a |
| unmodeled | z-up | 20 | 1.015x | 1.201x | 0.982x | 1.018x | 0.987x | 0.0772 | 0.0772 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.001x | 1.003x | 1.008x | 0.986x | 0.0772 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.003x | 0.988x | 1.023x | 0.972x | 0.0772 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.998x | 1.005x | 0.987x | 1.023x | 0.974x | 0.0772 | n/a | n/a |
| unmodeled | z-up | 30 | 1.028x | 1.370x | 0.974x | 1.047x | 0.982x | 0.1158 | 0.1158 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1158 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.001x | 1.002x | 1.005x | 1.012x | 0.981x | 0.1158 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.004x | 0.984x | 1.056x | 0.961x | 0.1158 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.009x | 0.983x | 1.055x | 0.963x | 0.1158 | n/a | n/a |
| unmodeled | z-up | 50 | 1.064x | 1.785x | 0.958x | 1.114x | 0.977x | 0.1930 | 0.1930 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1930 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.001x | 1.005x | 1.009x | 1.021x | 0.977x | 0.1930 | n/a | 98% |
| mekf | z-up | 50 | 1.002x | 1.006x | 0.974x | 1.125x | 0.953x | 0.1930 | n/a | 97% |
| mekf_ou3 | z-up | 50 | 0.997x | 1.018x | 0.973x | 1.128x | 0.947x | 0.1930 | n/a | 105% |
| unmodeled | z-up | 100 | 1.207x | 2.999x | 0.923x | 1.345x | 0.993x | 0.3860 | 0.3860 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3860 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.015x | 1.017x | 1.037x | 1.031x | 0.3860 | n/a | 99% |
| mekf | z-up | 100 | 1.005x | 1.016x | 0.956x | 1.371x | 1.038x | 0.3860 | n/a | 98% |
| mekf_ou3 | z-up | 100 | 0.996x | 1.054x | 0.953x | 1.378x | 0.952x | 0.3860 | n/a | 102% |
| unmodeled | z-up | 200 | 1.638x | 5.610x | 0.886x | 1.818x | 1.164x | 0.7720 | 0.7720 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7720 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.021x | 1.065x | 1.030x | 1.079x | 1.175x | 0.7720 | n/a | 97% |
| mekf | z-up | 200 | 1.013x | 1.052x | 0.947x | 1.868x | 1.764x | 0.7720 | n/a | 98% |
| mekf_ou3 | z-up | 200 | 1.005x | 1.183x | 0.948x | 1.892x | 1.228x | 0.7720 | n/a | 99% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2125 m, max roll/pitch 0.2554 deg, yaw 0.5819 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.001x | 1.103x | 1.005x | 1.013x | 0.999x | 0.0258 | 0.0258 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0258 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.002x | 0.998x | 1.003x | 1.000x | 0.0258 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.018x | 1.014x | 1.034x | 0.999x | 0.0258 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.001x | 1.022x | 1.073x | 1.142x | 0.991x | 0.0258 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.004x | 1.218x | 1.010x | 1.023x | 0.998x | 0.0516 | 0.0516 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0516 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 0.999x | 1.003x | 0.995x | 1.005x | 1.001x | 0.0516 | n/a | n/a |
| mekf | x-forward | 20 | 1.001x | 1.027x | 1.028x | 1.068x | 0.997x | 0.0516 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.003x | 1.034x | 1.154x | 1.297x | 0.987x | 0.0516 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.009x | 1.354x | 1.010x | 1.025x | 0.998x | 0.0773 | 0.0773 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0773 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 0.999x | 1.005x | 0.993x | 1.007x | 1.002x | 0.0773 | n/a | n/a |
| mekf | x-forward | 30 | 1.001x | 1.042x | 1.038x | 1.090x | 0.997x | 0.0773 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.006x | 1.053x | 1.236x | 1.447x | 0.989x | 0.0773 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.025x | 1.648x | 1.012x | 1.042x | 0.997x | 0.1289 | 0.1289 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1289 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 0.999x | 1.007x | 0.989x | 1.010x | 1.004x | 0.1289 | n/a | n/a |
| mekf | x-forward | 50 | 1.004x | 1.055x | 1.066x | 1.171x | 0.995x | 0.1289 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.015x | 1.078x | 1.426x | 1.798x | 1.003x | 0.1289 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.092x | 2.500x | 1.005x | 1.087x | 1.009x | 0.2578 | 0.2578 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2578 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 0.999x | 1.013x | 0.978x | 1.017x | 1.010x | 0.2578 | n/a | 101% |
| mekf | x-forward | 100 | 1.012x | 1.150x | 1.151x | 1.354x | 1.003x | 0.2578 | n/a | 87% |
| mekf_ou3 | x-forward | 100 | 1.045x | 1.215x | 1.995x | 2.796x | 1.264x | 0.2578 | n/a | 51% |
| unmodeled | x-forward | 200 | 1.318x | 4.269x | 0.940x | 1.187x | 1.115x | 0.5155 | 0.5155 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5155 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.016x | 1.067x | 0.974x | 1.103x | 1.013x | 0.5155 | n/a | 95% |
| mekf | x-forward | 200 | 1.034x | 1.388x | 1.303x | 1.717x | 1.075x | 0.5155 | n/a | 89% |
| mekf_ou3 | x-forward | 200 | 1.176x | 1.543x | 2.482x | 3.462x | 11.727x | 0.5155 | n/a | 45% |
| unmodeled | y-starboard | 10 | 1.004x | 1.116x | 0.995x | 1.009x | 1.008x | 0.0298 | 0.0298 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0298 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.001x | 0.999x | 1.001x | 1.000x | 0.0298 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.006x | 1.001x | 1.009x | 0.997x | 0.0298 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 1.000x | 1.005x | 1.010x | 1.028x | 0.946x | 0.0298 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.011x | 1.256x | 0.991x | 1.018x | 1.018x | 0.0597 | 0.0597 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0597 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.001x | 1.001x | 0.999x | 1.002x | 1.000x | 0.0597 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.009x | 1.002x | 1.019x | 0.994x | 0.0597 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 1.001x | 1.010x | 1.021x | 1.054x | 0.918x | 0.0597 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.019x | 1.421x | 0.985x | 1.027x | 1.028x | 0.0895 | 0.0895 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0895 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.001x | 1.002x | 0.999x | 1.002x | 1.000x | 0.0895 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.005x | 1.002x | 1.025x | 0.990x | 0.0895 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 1.002x | 1.008x | 1.032x | 1.080x | 0.921x | 0.0895 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.043x | 1.793x | 0.972x | 1.045x | 1.054x | 0.1492 | 0.1492 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1492 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.003x | 1.003x | 0.999x | 1.002x | 1.000x | 0.1492 | n/a | n/a |
| mekf | y-starboard | 50 | 1.001x | 1.003x | 1.000x | 1.030x | 0.985x | 0.1492 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 1.008x | 1.009x | 1.055x | 1.132x | 1.027x | 0.1492 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.137x | 2.813x | 0.935x | 1.112x | 1.134x | 0.2984 | 0.2984 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.2984 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.007x | 1.008x | 0.993x | 0.999x | 1.005x | 0.2984 | n/a | 95% |
| mekf | y-starboard | 100 | 1.006x | 1.017x | 0.989x | 1.035x | 0.978x | 0.2984 | n/a | 95% |
| mekf_ou3 | y-starboard | 100 | 1.033x | 1.041x | 1.217x | 2.183x | 1.872x | 0.2984 | n/a | 76% |
| unmodeled | y-starboard | 200 | 1.425x | 4.976x | 0.989x | 1.209x | 1.231x | 0.5968 | 0.5968 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.5968 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.032x | 1.037x | 0.989x | 1.047x | 1.015x | 0.5968 | n/a | 93% |
| mekf | y-starboard | 200 | 1.031x | 1.038x | 0.930x | 0.975x | 1.008x | 0.5968 | n/a | 93% |
| mekf_ou3 | y-starboard | 200 | 1.138x | 1.183x | 2.536x | 5.053x | 8.385x | 0.5968 | n/a | 67% |
| unmodeled | z-up | 10 | 1.006x | 1.072x | 0.991x | 1.010x | 0.993x | 0.0386 | 0.0386 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0386 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.001x | 1.001x | 1.002x | 0.996x | 0.0386 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.002x | 0.993x | 1.006x | 0.988x | 0.0386 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.999x | 1.003x | 0.993x | 1.009x | 0.990x | 0.0386 | n/a | n/a |
| unmodeled | z-up | 20 | 1.015x | 1.201x | 0.983x | 1.022x | 0.986x | 0.0772 | 0.0772 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0772 | 0.0000 | n/a |
| wrapper | z-up | 20 | 1.000x | 1.001x | 1.002x | 1.004x | 0.992x | 0.0772 | n/a | n/a |
| mekf | z-up | 20 | 1.001x | 1.003x | 0.986x | 1.024x | 0.977x | 0.0772 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.998x | 1.005x | 0.987x | 1.025x | 0.981x | 0.0772 | n/a | n/a |
| unmodeled | z-up | 30 | 1.028x | 1.370x | 0.977x | 1.051x | 0.980x | 0.1158 | 0.1158 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1158 | 0.0000 | n/a |
| wrapper | z-up | 30 | 1.001x | 1.002x | 1.003x | 1.006x | 0.989x | 0.1158 | n/a | n/a |
| mekf | z-up | 30 | 1.001x | 1.004x | 0.981x | 1.054x | 0.967x | 0.1158 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.998x | 1.009x | 0.982x | 1.056x | 0.972x | 0.1158 | n/a | n/a |
| unmodeled | z-up | 50 | 1.064x | 1.785x | 0.962x | 1.115x | 0.972x | 0.1930 | 0.1930 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1930 | 0.0000 | 100% |
| wrapper | z-up | 50 | 1.001x | 1.005x | 1.005x | 1.010x | 0.983x | 0.1930 | n/a | 99% |
| mekf | z-up | 50 | 1.002x | 1.006x | 0.970x | 1.119x | 0.952x | 0.1930 | n/a | 97% |
| mekf_ou3 | z-up | 50 | 0.997x | 1.018x | 0.970x | 1.125x | 0.958x | 0.1930 | n/a | 105% |
| unmodeled | z-up | 100 | 1.207x | 3.000x | 0.932x | 1.346x | 0.974x | 0.3860 | 0.3860 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3860 | 0.0000 | 100% |
| wrapper | z-up | 100 | 1.002x | 1.015x | 1.010x | 1.021x | 0.988x | 0.3860 | n/a | 99% |
| mekf | z-up | 100 | 1.005x | 1.016x | 0.947x | 1.355x | 0.960x | 0.3860 | n/a | 98% |
| mekf_ou3 | z-up | 100 | 0.996x | 1.054x | 0.949x | 1.369x | 0.945x | 0.3860 | n/a | 102% |
| unmodeled | z-up | 200 | 1.638x | 5.612x | 0.904x | 1.819x | 1.077x | 0.7720 | 0.7720 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.7720 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.021x | 1.065x | 1.022x | 1.063x | 1.046x | 0.7720 | n/a | 97% |
| mekf | z-up | 200 | 1.012x | 1.052x | 0.936x | 1.838x | 1.286x | 0.7720 | n/a | 98% |
| mekf_ou3 | z-up | 200 | 1.005x | 1.183x | 0.941x | 1.869x | 1.053x | 0.7720 | n/a | 99% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.129x | 1.033x | 1.031x |
| CARMA3 | 2 | 1.030x | 1.026x | 1.018x |
| CARMA3 | 5 | 1.007x | 1.014x | 0.999x |
| CARMA3 | 10 | 1.003x | 1.013x | 0.982x |
| CARMA3 | 15 | 1.002x | 1.017x | 1.031x |
| CARMA3 | 25 | 1.021x | 1.025x | 1.148x |
| CARMA3 | 50 | 1.114x | 1.032x | 1.007x |
| CARMA3 | 100 | 1.106x | 1.027x | 1.004x |
| TFG-CARMA3 | 1 | 1.129x | 1.036x | 1.025x |
| TFG-CARMA3 | 2 | 1.030x | 1.027x | 1.017x |
| TFG-CARMA3 | 5 | 1.007x | 1.013x | 1.003x |
| TFG-CARMA3 | 10 | 1.003x | 1.009x | 0.987x |
| TFG-CARMA3 | 15 | 1.002x | 1.010x | 0.988x |
| TFG-CARMA3 | 25 | 1.021x | 1.017x | 1.026x |
| TFG-CARMA3 | 50 | 1.114x | 1.031x | 0.999x |
| TFG-CARMA3 | 100 | 1.106x | 1.027x | 0.995x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
