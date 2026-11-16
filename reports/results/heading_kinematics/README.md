# Heading and horizontal-wave attribution

[Analysis, limitations and reproduction](../../../docs/heading-kinematics.md).

- `measurements.json` / `.csv`: 48 low-wave cells, eight families, two directions, three physical histories.
- `medium.json` / `.csv`: 12 medium-wave confirmation cells for OU-III and shared QMEKF.
- SVGs and `*-tail.csv`: final 60 s, uniformly sampled at 20 Hz; numerical scores use the full 900 s / 200 Hz window.

- `runtime.json` / `.csv`: the implemented turn model, all 128 default-profile cells, with unchanged gate failures retained.
- `runtime-controls.json`: 20 zero-speed controls of that runtime.
- `runtime-native-ci.json`: the separate native-target matrix, both directions, with failed profiles and complete CI artifact identities preserved.
- `runtime-small-ou3.svg` / `runtime-*-small-ou3-tail.csv`: pre-fix versus runtime, same sensor stream; the slow bias offset remains visible.

The original attribution bundles and runtime measurements have separate roles.
The runtime reduces Y RMS in all 128 cells and the measured bias harmonic in
all 20 paired cells, but only 15/128 cells pass every existing quality gate.
These files do not replace the locked acceptance baseline.
