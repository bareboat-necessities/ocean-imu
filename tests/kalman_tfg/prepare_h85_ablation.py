#!/usr/bin/env python3
"""Prepare one targeted TFG Hs=8.5 ablation in a GitHub Actions workspace.

The committed production source is the corrected implementation.  This helper
edits the checked-out workspace only, so every arm starts from the same branch
and changes exactly one correction family.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
KALMAN = ROOT / "src/kalman_tfg/Kalman3D_Wave_TFG.h"
FUSION = ROOT / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {n}")
    return text.replace(old, new, 1)


def patch_joseph_only(text: str) -> str:
    old = """        MatrixNX& Jreset = scratch_c_;\n        build_reset_jacobian(correction, Jreset);\n        inject(Tangent(-correction));\n        transport_reset_(Jreset, Pj);\n"""
    new = """        // Ablation: old Joseph-only covariance handling.  The nominal\n        // state is still retracted, but P is intentionally left in the\n        // pre-reset tangent coordinates.\n        inject(Tangent(-correction));\n        P_ = Pj;\n        symmetrize_in_place_(P_);\n"""
    return replace_once(text, old, new, "Joseph-only reset ablation")


def patch_old_mag_lock(text: str) -> str:
    old = """                const float psi = std::atan2(b_world.y(), b_world.x());\n                mekf_.apply_world_yaw_gauge(-psi);\n"""
    new = """                const float psi = std::atan2(b_world.y(), b_world.x());\n                // Ablation: historical first-lock behaviour.  Rotate only R;\n                // do not transform world columns or covariance.\n                const Eigen::Matrix3f Rz =\n                    Eigen::AngleAxisf(-psi, Vector3f::UnitZ()).toRotationMatrix();\n                mekf_.state().R = (Rz * mekf_.state().R).eval();\n"""
    return replace_once(text, old, new, "old magnetic lock ablation")


def patch_old_bias_qd(text: str) -> str:
    # A random-walk bias has constant mean and perturbation magnitude. Keep
    # the existing endpoint rotation of the Level-2 bias coordinates.
    text = replace_once(text, "X_.B.col(1) *= std::exp(-Ts / tau_ba_);",
                        "// Historical random walk: constant bias.",
                        "standalone bias mean")
    text = replace_once(text,
                        "s.alpha_ba = with_accel_bias && acc_bias_updates_enabled_ ? std::exp(-Ts/tau_ba_) : T(1);",
                        "s.alpha_ba = T(1);", "propagated bias mean and transition")
    text = replace_once(text,
                        "const T qint = T(0.5) * tau_ba_ * (-std::expm1(T(-2)*Ts/tau_ba_));",
                        "const T qint = Ts;", "bias random-walk covariance")

    # Historical constant-rotation h^2/2 and h^3/3 gyro-bias integrals.
    historical_noise = """        const Matrix3 Wb = X_.R * Q_bg_ * X_.R.transpose();
        const Matrix3 E = two_frame_bias ? step.R_end : Matrix3::Identity();
        Qtt = (X_.R * (gyro_noise_scale * Q_gyro_) * X_.R.transpose()) * Ts
            + Wb * (Ts * Ts * Ts / T(3));
        Qtb = -(X_.R * Q_bg_ * E.transpose()) * (Ts * Ts / T(2));
        Qbb = (E * Q_bg_ * E.transpose()) * Ts;"""
    return replace_once(text, "        gyro_noise_(Ts, step, gyro_noise_scale, Qtt, Qtb, Qbb);",
                        historical_noise, "historical gyro-bias process covariance")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: prepare_h85_ablation.py corrected|joseph_only|old_mag_lock|old_bias_qd",
              file=sys.stderr)
        return 2

    variant = sys.argv[1]
    allowed = {"corrected", "joseph_only", "old_mag_lock", "old_bias_qd"}
    if variant not in allowed:
        print(f"unknown variant {variant!r}", file=sys.stderr)
        return 2

    kalman = KALMAN.read_text()
    fusion = FUSION.read_text()

    if variant == "joseph_only":
        kalman = patch_joseph_only(kalman)
    elif variant == "old_mag_lock":
        fusion = patch_old_mag_lock(fusion)
    elif variant == "old_bias_qd":
        kalman = patch_old_bias_qd(kalman)

    KALMAN.write_text(kalman)
    FUSION.write_text(fusion)
    print(f"prepared TFG H8.5 ablation variant: {variant}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
