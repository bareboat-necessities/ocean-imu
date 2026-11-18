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

CG baseline: 3-D displacement 0.2286 m, max roll/pitch 0.1804 deg, yaw 0.7701 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.000x | 1.065x | 0.941x | 0.997x | 1.010x | 0.0166 | 0.0166 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0166 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 1.002x | 1.004x | 0.981x | 0.0166 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.004x | 1.004x | 1.008x | 0.981x | 0.0166 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.008x | 1.076x | 1.120x | 1.008x | 0.0166 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.001x | 1.159x | 0.889x | 0.993x | 1.020x | 0.0331 | 0.0331 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.001x | 0.0331 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.002x | 1.004x | 1.008x | 0.965x | 0.0331 | n/a | n/a |
| mekf | x-forward | 20 | 1.000x | 1.010x | 1.009x | 1.017x | 0.963x | 0.0331 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.020x | 1.158x | 1.256x | 1.018x | 0.0331 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.004x | 1.273x | 0.847x | 0.990x | 1.030x | 0.0497 | 0.0497 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0497 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.003x | 1.005x | 1.012x | 0.950x | 0.0497 | n/a | n/a |
| mekf | x-forward | 30 | 1.000x | 1.019x | 1.014x | 1.026x | 0.946x | 0.0497 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.006x | 1.034x | 1.248x | 1.406x | 1.029x | 0.0497 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.015x | 1.541x | 0.793x | 0.983x | 1.051x | 0.0828 | 0.0828 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0828 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 1.000x | 1.006x | 1.011x | 1.023x | 0.930x | 0.0828 | n/a | n/a |
| mekf | x-forward | 50 | 1.001x | 1.039x | 1.025x | 1.044x | 0.922x | 0.0828 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.013x | 1.069x | 1.443x | 1.728x | 1.056x | 0.0828 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.068x | 2.281x | 0.891x | 1.107x | 1.107x | 0.1657 | 0.1657 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.001x | 0.1657 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 1.002x | 1.016x | 1.027x | 1.057x | 0.954x | 0.1657 | n/a | 98% |
| mekf | x-forward | 100 | 1.005x | 1.075x | 1.051x | 1.074x | 0.938x | 0.1657 | n/a | 93% |
| mekf_ou3 | x-forward | 100 | 1.052x | 1.188x | 2.312x | 2.989x | 2.040x | 0.1657 | n/a | 23% |
| unmodeled | x-forward | 200 | 1.264x | 3.801x | 1.650x | 3.182x | 1.204x | 0.3313 | 0.3313 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.3313 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.009x | 1.064x | 1.088x | 1.137x | 1.133x | 0.3313 | n/a | 97% |
| mekf | x-forward | 200 | 1.020x | 1.066x | 1.081x | 1.158x | 1.383x | 0.3313 | n/a | 92% |
| mekf_ou3 | x-forward | 200 | 1.225x | 1.582x | 4.775x | 7.180x | 5.064x | 0.3313 | n/a | 15% |
| unmodeled | y-starboard | 10 | 1.003x | 1.009x | 0.995x | 1.001x | 1.044x | 0.0088 | 0.0088 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0088 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.003x | 1.005x | 0.988x | 0.0088 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.000x | 1.003x | 1.006x | 0.985x | 0.0088 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 0.998x | 1.010x | 1.016x | 0.957x | 0.0088 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.008x | 1.026x | 0.989x | 1.002x | 1.093x | 0.0175 | 0.0175 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0175 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.005x | 1.009x | 0.974x | 0.0175 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.000x | 1.006x | 1.011x | 0.974x | 0.0175 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.997x | 0.997x | 1.019x | 1.033x | 0.955x | 0.0175 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.015x | 1.070x | 0.984x | 1.159x | 1.143x | 0.0263 | 0.0263 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0263 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.007x | 1.012x | 0.964x | 0.0263 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.000x | 1.008x | 1.016x | 0.964x | 0.0263 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.996x | 0.997x | 1.028x | 1.049x | 0.992x | 0.0263 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.032x | 1.211x | 1.108x | 1.577x | 1.250x | 0.0438 | 0.0438 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0438 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 0.999x | 1.000x | 1.009x | 1.018x | 0.948x | 0.0438 | n/a | n/a |
| mekf | y-starboard | 50 | 1.000x | 1.000x | 1.011x | 1.024x | 0.948x | 0.0438 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.996x | 1.005x | 1.051x | 1.216x | 1.177x | 0.0438 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.097x | 1.750x | 1.663x | 2.803x | 1.538x | 0.0877 | 0.0877 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0877 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.000x | 1.002x | 1.008x | 1.023x | 0.957x | 0.0877 | n/a | 100% |
| mekf | y-starboard | 100 | 1.001x | 1.003x | 1.009x | 1.034x | 0.946x | 0.0877 | n/a | 99% |
| mekf_ou3 | y-starboard | 100 | 1.015x | 1.133x | 1.695x | 2.458x | 2.650x | 0.0877 | n/a | 84% |
| unmodeled | y-starboard | 200 | 1.305x | 3.118x | 2.842x | 5.203x | 2.158x | 0.1754 | 0.1754 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.1754 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.003x | 1.032x | 0.988x | 1.013x | 1.047x | 0.1754 | n/a | 99% |
| mekf | y-starboard | 200 | 1.007x | 1.074x | 0.981x | 1.022x | 1.396x | 0.1754 | n/a | 98% |
| mekf_ou3 | y-starboard | 200 | 1.201x | 1.669x | 6.129x | 8.001x | 8.796x | 0.1754 | n/a | 34% |
| unmodeled | z-up | 10 | 1.011x | 1.041x | 1.000x | 1.002x | 1.001x | 0.0185 | 0.0185 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0185 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.004x | 0.976x | 0.0185 | n/a | n/a |
| mekf | z-up | 10 | 1.000x | 1.000x | 1.002x | 1.005x | 0.976x | 0.0185 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.998x | 0.999x | 1.001x | 1.004x | 0.986x | 0.0185 | n/a | n/a |
| unmodeled | z-up | 20 | 1.023x | 1.121x | 1.000x | 1.004x | 1.001x | 0.0370 | 0.0370 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0370 | 0.0000 | n/a |
| wrapper | z-up | 20 | 0.999x | 1.000x | 1.005x | 1.007x | 0.947x | 0.0370 | n/a | n/a |
| mekf | z-up | 20 | 1.000x | 1.000x | 1.005x | 1.012x | 0.947x | 0.0370 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.997x | 0.997x | 1.003x | 1.009x | 0.968x | 0.0370 | n/a | n/a |
| unmodeled | z-up | 30 | 1.038x | 1.232x | 0.999x | 1.006x | 0.999x | 0.0556 | 0.0556 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0556 | 0.0000 | n/a |
| wrapper | z-up | 30 | 0.999x | 1.000x | 1.008x | 1.012x | 0.918x | 0.0556 | n/a | n/a |
| mekf | z-up | 30 | 0.999x | 1.000x | 1.008x | 1.018x | 0.917x | 0.0556 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.995x | 0.996x | 1.004x | 1.014x | 0.949x | 0.0556 | n/a | n/a |
| unmodeled | z-up | 50 | 1.075x | 1.519x | 0.999x | 1.010x | 1.001x | 0.0926 | 0.0926 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.0926 | 0.0000 | 100% |
| wrapper | z-up | 50 | 0.998x | 1.001x | 1.015x | 1.022x | 0.868x | 0.0926 | n/a | 102% |
| mekf | z-up | 50 | 0.999x | 1.001x | 1.015x | 1.033x | 0.866x | 0.0926 | n/a | 102% |
| mekf_ou3 | z-up | 50 | 0.992x | 0.995x | 1.008x | 1.024x | 0.913x | 0.0926 | n/a | 111% |
| unmodeled | z-up | 100 | 1.201x | 2.414x | 0.999x | 1.021x | 1.002x | 0.1852 | 0.1852 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.001x | 1.001x | 0.1852 | 0.0000 | 100% |
| wrapper | z-up | 100 | 0.998x | 1.004x | 1.028x | 1.046x | 0.896x | 0.1852 | n/a | 101% |
| mekf | z-up | 100 | 0.998x | 1.005x | 1.028x | 1.066x | 0.883x | 0.1852 | n/a | 101% |
| mekf_ou3 | z-up | 100 | 0.985x | 0.997x | 1.018x | 1.052x | 0.837x | 0.1852 | n/a | 107% |
| unmodeled | z-up | 200 | 1.549x | 4.416x | 1.003x | 1.047x | 1.007x | 0.3704 | 0.3704 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.001x | 0.3704 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.001x | 1.033x | 1.015x | 1.041x | 1.183x | 0.3704 | n/a | 100% |
| mekf | z-up | 200 | 1.015x | 1.214x | 1.113x | 2.429x | 5.068x | 0.3704 | n/a | 97% |
| mekf_ou3 | z-up | 200 | 0.977x | 1.024x | 1.047x | 1.092x | 1.210x | 0.3704 | n/a | 104% |

## TFG-CARMA3

CG baseline: 3-D displacement 0.2252 m, max roll/pitch 0.1788 deg, yaw 0.7628 deg.

| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea | Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| unmodeled | x-forward | 10 | 1.000x | 1.065x | 0.940x | 0.996x | 1.010x | 0.0166 | 0.0166 | n/a |
| exact | x-forward | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0166 | 0.0000 | n/a |
| wrapper | x-forward | 10 | 1.000x | 1.001x | 1.000x | 1.001x | 0.990x | 0.0166 | n/a | n/a |
| mekf | x-forward | 10 | 1.000x | 1.004x | 1.003x | 1.009x | 0.989x | 0.0166 | n/a | n/a |
| mekf_ou3 | x-forward | 10 | 1.002x | 1.008x | 1.080x | 1.122x | 0.997x | 0.0166 | n/a | n/a |
| unmodeled | x-forward | 20 | 1.001x | 1.159x | 0.887x | 0.992x | 1.019x | 0.0331 | 0.0331 | n/a |
| exact | x-forward | 20 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.0331 | 0.0000 | n/a |
| wrapper | x-forward | 20 | 1.000x | 1.002x | 1.001x | 1.006x | 0.980x | 0.0331 | n/a | n/a |
| mekf | x-forward | 20 | 1.000x | 1.010x | 1.006x | 1.018x | 0.979x | 0.0331 | n/a | n/a |
| mekf_ou3 | x-forward | 20 | 1.004x | 1.020x | 1.167x | 1.262x | 0.995x | 0.0331 | n/a | n/a |
| unmodeled | x-forward | 30 | 1.005x | 1.273x | 0.845x | 0.988x | 1.030x | 0.0497 | 0.0497 | n/a |
| exact | x-forward | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0497 | 0.0000 | n/a |
| wrapper | x-forward | 30 | 1.000x | 1.003x | 1.001x | 1.003x | 0.972x | 0.0497 | n/a | n/a |
| mekf | x-forward | 30 | 1.000x | 1.019x | 1.009x | 1.027x | 0.969x | 0.0497 | n/a | n/a |
| mekf_ou3 | x-forward | 30 | 1.006x | 1.035x | 1.261x | 1.415x | 0.996x | 0.0497 | n/a | n/a |
| unmodeled | x-forward | 50 | 1.016x | 1.541x | 0.795x | 0.981x | 1.051x | 0.0828 | 0.0828 | n/a |
| exact | x-forward | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0828 | 0.0000 | n/a |
| wrapper | x-forward | 50 | 1.000x | 1.006x | 1.002x | 1.006x | 0.955x | 0.0828 | n/a | n/a |
| mekf | x-forward | 50 | 1.001x | 1.039x | 1.016x | 1.045x | 0.950x | 0.0828 | n/a | n/a |
| mekf_ou3 | x-forward | 50 | 1.012x | 1.069x | 1.462x | 1.740x | 0.997x | 0.0828 | n/a | n/a |
| unmodeled | x-forward | 100 | 1.071x | 2.281x | 0.912x | 1.285x | 1.108x | 0.1657 | 0.1657 | n/a |
| exact | x-forward | 100 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.1657 | 0.0000 | 100% |
| wrapper | x-forward | 100 | 1.001x | 1.016x | 1.005x | 1.016x | 0.929x | 0.1657 | n/a | 98% |
| mekf | x-forward | 100 | 1.004x | 1.075x | 1.029x | 1.072x | 0.914x | 0.1657 | n/a | 95% |
| mekf_ou3 | x-forward | 100 | 1.040x | 1.131x | 1.810x | 2.900x | 2.006x | 0.1657 | n/a | 44% |
| unmodeled | x-forward | 200 | 1.272x | 3.800x | 1.711x | 3.534x | 1.208x | 0.3313 | 0.3313 | n/a |
| exact | x-forward | 200 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.3313 | 0.0000 | 100% |
| wrapper | x-forward | 200 | 1.006x | 1.064x | 1.030x | 1.054x | 0.842x | 0.3313 | n/a | 98% |
| mekf | x-forward | 200 | 1.014x | 1.065x | 1.028x | 1.064x | 0.934x | 0.3313 | n/a | 95% |
| mekf_ou3 | x-forward | 200 | 1.838x | 4.675x | 5.111x | 9.822x | 54.788x | 0.3313 | n/a | -208% |
| unmodeled | y-starboard | 10 | 1.003x | 1.009x | 0.995x | 1.001x | 1.047x | 0.0088 | 0.0088 | n/a |
| exact | y-starboard | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0088 | 0.0000 | n/a |
| wrapper | y-starboard | 10 | 1.000x | 1.000x | 1.001x | 1.006x | 0.993x | 0.0088 | n/a | n/a |
| mekf | y-starboard | 10 | 1.000x | 1.000x | 1.001x | 1.003x | 0.993x | 0.0088 | n/a | n/a |
| mekf_ou3 | y-starboard | 10 | 0.998x | 0.998x | 1.008x | 1.015x | 0.952x | 0.0088 | n/a | n/a |
| unmodeled | y-starboard | 20 | 1.008x | 1.026x | 0.989x | 1.002x | 1.097x | 0.0175 | 0.0175 | n/a |
| exact | y-starboard | 20 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0175 | 0.0000 | n/a |
| wrapper | y-starboard | 20 | 1.000x | 1.000x | 1.002x | 1.003x | 0.986x | 0.0175 | n/a | n/a |
| mekf | y-starboard | 20 | 1.000x | 1.000x | 1.002x | 1.006x | 0.986x | 0.0175 | n/a | n/a |
| mekf_ou3 | y-starboard | 20 | 0.997x | 0.997x | 1.016x | 1.029x | 0.924x | 0.0175 | n/a | n/a |
| unmodeled | y-starboard | 30 | 1.015x | 1.071x | 0.984x | 1.156x | 1.150x | 0.0263 | 0.0263 | n/a |
| exact | y-starboard | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0263 | 0.0000 | n/a |
| wrapper | y-starboard | 30 | 1.000x | 1.000x | 1.002x | 1.004x | 0.980x | 0.0263 | n/a | n/a |
| mekf | y-starboard | 30 | 1.000x | 1.000x | 1.002x | 1.008x | 0.979x | 0.0263 | n/a | n/a |
| mekf_ou3 | y-starboard | 30 | 0.996x | 0.996x | 1.024x | 1.044x | 0.925x | 0.0263 | n/a | n/a |
| unmodeled | y-starboard | 50 | 1.032x | 1.211x | 1.121x | 1.655x | 1.261x | 0.0438 | 0.0438 | n/a |
| exact | y-starboard | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0438 | 0.0000 | n/a |
| wrapper | y-starboard | 50 | 1.000x | 1.000x | 1.002x | 1.005x | 0.968x | 0.0438 | n/a | n/a |
| mekf | y-starboard | 50 | 0.999x | 1.000x | 1.003x | 1.012x | 0.967x | 0.0438 | n/a | n/a |
| mekf_ou3 | y-starboard | 50 | 0.995x | 0.997x | 1.033x | 1.075x | 1.051x | 0.0438 | n/a | n/a |
| unmodeled | y-starboard | 100 | 1.099x | 1.751x | 1.701x | 2.993x | 1.564x | 0.0877 | 0.0877 | n/a |
| exact | y-starboard | 100 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0877 | 0.0000 | 100% |
| wrapper | y-starboard | 100 | 1.000x | 1.001x | 0.996x | 1.002x | 0.946x | 0.0877 | n/a | 100% |
| mekf | y-starboard | 100 | 1.000x | 1.003x | 0.997x | 1.014x | 0.944x | 0.0877 | n/a | 100% |
| mekf_ou3 | y-starboard | 100 | 1.011x | 1.063x | 1.420x | 2.198x | 2.463x | 0.0877 | n/a | 89% |
| unmodeled | y-starboard | 200 | 1.311x | 3.119x | 2.935x | 5.660x | 2.217x | 0.1754 | 0.1754 | n/a |
| exact | y-starboard | 200 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.1754 | 0.0000 | 100% |
| wrapper | y-starboard | 200 | 1.003x | 1.032x | 0.980x | 1.002x | 0.958x | 0.1754 | n/a | 99% |
| mekf | y-starboard | 200 | 1.002x | 1.074x | 0.966x | 0.992x | 0.944x | 0.1754 | n/a | 99% |
| mekf_ou3 | y-starboard | 200 | 3.188x | 8.083x | 5.467x | 9.387x | 67.590x | 0.1754 | n/a | -604% |
| unmodeled | z-up | 10 | 1.010x | 1.041x | 1.000x | 1.003x | 1.001x | 0.0185 | 0.0185 | n/a |
| exact | z-up | 10 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0185 | 0.0000 | n/a |
| wrapper | z-up | 10 | 1.000x | 1.000x | 1.001x | 1.002x | 0.987x | 0.0185 | n/a | n/a |
| mekf | z-up | 10 | 0.999x | 1.000x | 1.001x | 1.004x | 0.988x | 0.0185 | n/a | n/a |
| mekf_ou3 | z-up | 10 | 0.998x | 0.998x | 1.000x | 1.004x | 0.993x | 0.0185 | n/a | n/a |
| unmodeled | z-up | 20 | 1.024x | 1.121x | 0.999x | 1.005x | 1.001x | 0.0370 | 0.0370 | n/a |
| exact | z-up | 20 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0370 | 0.0000 | n/a |
| wrapper | z-up | 20 | 0.999x | 1.000x | 1.002x | 1.003x | 0.972x | 0.0370 | n/a | n/a |
| mekf | z-up | 20 | 0.999x | 1.000x | 1.001x | 1.009x | 0.973x | 0.0370 | n/a | n/a |
| mekf_ou3 | z-up | 20 | 0.996x | 0.997x | 1.000x | 1.008x | 0.985x | 0.0370 | n/a | n/a |
| unmodeled | z-up | 30 | 1.039x | 1.232x | 0.998x | 1.007x | 1.001x | 0.0556 | 0.0556 | n/a |
| exact | z-up | 30 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0556 | 0.0000 | n/a |
| wrapper | z-up | 30 | 0.999x | 1.000x | 1.003x | 1.004x | 0.956x | 0.0556 | n/a | n/a |
| mekf | z-up | 30 | 0.999x | 1.000x | 1.002x | 1.012x | 0.957x | 0.0556 | n/a | n/a |
| mekf_ou3 | z-up | 30 | 0.995x | 0.996x | 1.000x | 1.011x | 0.976x | 0.0556 | n/a | n/a |
| unmodeled | z-up | 50 | 1.077x | 1.518x | 0.997x | 1.012x | 1.001x | 0.0926 | 0.0926 | n/a |
| exact | z-up | 50 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.0926 | 0.0000 | 100% |
| wrapper | z-up | 50 | 0.998x | 1.001x | 1.005x | 1.007x | 0.921x | 0.0926 | n/a | 102% |
| mekf | z-up | 50 | 0.998x | 1.001x | 1.003x | 1.021x | 0.922x | 0.0926 | n/a | 102% |
| mekf_ou3 | z-up | 50 | 0.991x | 0.995x | 1.000x | 1.019x | 0.956x | 0.0926 | n/a | 111% |
| unmodeled | z-up | 100 | 1.204x | 2.411x | 0.993x | 1.019x | 1.003x | 0.1852 | 0.1852 | n/a |
| exact | z-up | 100 | 1.000x | 1.000x | 1.000x | 1.001x | 1.000x | 0.1852 | 0.0000 | 100% |
| wrapper | z-up | 100 | 0.997x | 1.004x | 1.009x | 1.013x | 0.868x | 0.1852 | n/a | 101% |
| mekf | z-up | 100 | 0.998x | 1.010x | 1.005x | 1.035x | 0.869x | 0.1852 | n/a | 101% |
| mekf_ou3 | z-up | 100 | 0.984x | 1.004x | 0.998x | 1.031x | 0.902x | 0.1852 | n/a | 108% |
| unmodeled | z-up | 200 | 1.559x | 4.412x | 0.993x | 1.036x | 1.006x | 0.3704 | 0.3704 | n/a |
| exact | z-up | 200 | 1.000x | 1.000x | 1.000x | 1.000x | 1.000x | 0.3704 | 0.0000 | 100% |
| wrapper | z-up | 200 | 1.001x | 1.033x | 1.007x | 1.030x | 0.994x | 0.3704 | n/a | 100% |
| mekf | z-up | 200 | 1.002x | 1.029x | 1.010x | 1.050x | 1.532x | 0.3704 | n/a | 100% |
| mekf_ou3 | z-up | 200 | 0.977x | 1.035x | 1.006x | 1.054x | 0.886x | 0.3704 | n/a | 104% |

## Derivative corner of the deployed correction

Wrapper arm, IMU 100 cm above the CG, pooled over the same seas.

| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |
| --- | ---: | ---: | ---: | ---: |
| CARMA3 | 1 | 1.024x | 1.018x | 0.997x |
| CARMA3 | 2 | 1.004x | 1.009x | 0.995x |
| CARMA3 | 5 | 0.999x | 1.007x | 0.966x |
| CARMA3 | 10 | 0.998x | 1.017x | 0.873x |
| CARMA3 | 15 | 0.998x | 1.028x | 0.896x |
| CARMA3 | 25 | 1.002x | 1.017x | 1.151x |
| CARMA3 | 50 | 1.016x | 1.018x | 1.023x |
| CARMA3 | 100 | 1.018x | 1.063x | 1.058x |
| TFG-CARMA3 | 1 | 1.023x | 1.014x | 0.998x |
| TFG-CARMA3 | 2 | 1.004x | 1.006x | 0.997x |
| TFG-CARMA3 | 5 | 0.998x | 1.003x | 0.982x |
| TFG-CARMA3 | 10 | 0.997x | 1.006x | 0.922x |
| TFG-CARMA3 | 15 | 0.997x | 1.009x | 0.868x |
| TFG-CARMA3 | 25 | 1.002x | 1.004x | 0.942x |
| TFG-CARMA3 | 50 | 1.019x | 1.020x | 1.021x |
| TFG-CARMA3 | 100 | 1.023x | 1.069x | 1.055x |

## Reading the tables

The installed force is deterministic, but it combines with each record's
existing estimator residual, so an arm can score slightly below its CG
baseline; that is not evidence that the installation helps.  The exact arm
bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets
the results.

Study matrix: 8 records, 2 filters, 3 axes, 6 offsets, 5 arms.
