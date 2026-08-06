# Vendored: Lie++ (Lie-plusplus)

Upstream: https://github.com/aau-cns/Lie-plusplus
Commit:   ec459c38445d986b5830162ae428b8b242ce5021
License:  Apache License 2.0 (see `LICENSE`)
Authors:  Alessandro Fornasier, Pieter van Goor et al., Control of Networked
          Systems (CNS), University of Klagenfurt.

## Why it is here

`tests/kalman_tfg/lie_group-test` uses `group::SEn3<double, 4>` as an
**independent oracle** for the world-vector block of the two-frame group
implemented in `src/lie/`. Checking against finite differences alone would
leave a shared-assumption blind spot; Lie++ is a separately derived
implementation by the authors of the equivariant filtering literature this
filter follows.

## Why it is not in `src/`

Lie++ is compiled **only** by native tests, never for Arduino:

- it includes `<Eigen/Dense>` directly, whereas firmware builds route through
  `ArduinoEigenDense.h` (see `src/ArduinoOceanImu.h`);
- it materializes dense `(3+n)x(3+n)` and `(3+3n)x(3+3n)` Eigen matrices, and
  `src/kalman_ou_iii/Kalman3D_Wave_OU_III.h:1707` documents why this codebase
  keeps fixed-size Eigen instantiations small (cc1plus OOM on ESP32/Windows).

The firmware path uses `src/lie/`, which is scalar-templated, Arduino-safe, and
written in the surrounding `ou_detail` style.

## Scope of the copy

Only the two headers the oracle needs:

- `include/groups/SO3.hpp`
- `include/groups/SEn3.hpp`  (depends on `SO3.hpp` only)

Upstream `SDB.hpp` is deliberately **not** vendored: its semi-direct bias
structure is hardcoded to `SE_2(3) + R^6`, so it cannot serve as an oracle for
this filter's four world vectors plus six bias states. That block is validated
by structural identities and finite differences instead.

Files are copied verbatim. Do not edit them; to update, re-copy from a newer
upstream commit and record it above.

## Convention agreement

Lie++'s conventions line up with `src/lie/TwoFrameGroup.h` directly, which is
what makes the comparison meaningful rather than an exercise in re-indexing:

| Operation | Lie++ `SEn3` | `TwoFrameGroup` |
|---|---|---|
| composition | `t_i = R1 t2_i + t1_i` | `X1 + R1 X2` |
| inverse | `-R^T t_i` | `-R^T X` |
| exponential | `t_i = J_l(w) u_i` | `E_X = J_l(phi) rho` |
| adjoint | `[t_i]x R` | `[x_i]x R` |

Tangent layouts differ: Lie++ packs `[w, t_0, t_1, ...]` contiguously, while
this filter interleaves the gyro bias at offset 3 to preserve OU-III's state
ordering. The world part is contiguous in both, so the mapping is a pair of
segment copies.

## Citation

A. Fornasier et al., "Equivariant Symmetries for Inertial Navigation Systems,"
arXiv:2309.03765, 2023.
