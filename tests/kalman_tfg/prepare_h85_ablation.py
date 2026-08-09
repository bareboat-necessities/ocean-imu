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


def replace_between(text: str, start: str, end: str, new: str, label: str) -> str:
    i = text.find(start)
    if i < 0:
        raise RuntimeError(f"{label}: start marker not found")
    j = text.find(end, i)
    if j < 0:
        raise RuntimeError(f"{label}: end marker not found")
    j += len(end)
    return text[:i] + new + text[j:]


def patch_joseph_only(text: str) -> str:
    old = """        MatrixNX Jreset;\n        build_reset_jacobian(correction, Jreset);\n        inject(Tangent(-correction));\n        P_ = (Jreset * Pj * Jreset.transpose()).eval();\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""
    new = """        // Ablation: old Joseph-only covariance handling.  The nominal\n        // state is still retracted, but P is intentionally left in the\n        // pre-reset tangent coordinates.\n        inject(Tangent(-correction));\n        P_ = Pj;\n        P_ = T(0.5) * (P_ + P_.transpose()).eval();\n"""
    return replace_once(text, old, new, "Joseph-only reset ablation")


def patch_old_mag_lock(text: str) -> str:
    old = """                const float psi = std::atan2(b_world.y(), b_world.x());\n                mekf_.apply_world_yaw_gauge(-psi);\n"""
    new = """                const float psi = std::atan2(b_world.y(), b_world.x());\n                // Ablation: historical first-lock behaviour.  Rotate only R;\n                // do not transform world columns or covariance.\n                const Eigen::Matrix3f Rz =\n                    Eigen::AngleAxisf(-psi, Vector3f::UnitZ()).toRotationMatrix();\n                mekf_.state().R = (Rz * mekf_.state().R).eval();\n"""
    return replace_once(text, old, new, "old magnetic lock ablation")


def patch_old_bias_qd(text: str) -> str:
    # 1) Residual accelerometer bias mean: OU -> historical random walk.
    old_mean = """\n        if constexpr (with_accel_bias) {\n            if (acc_bias_updates_enabled_) {\n                X_.B.col(1) *= std::exp(-Ts / tau_ba_);\n            }\n        }\n"""
    new_mean = """\n        // Ablation: historical accelerometer-bias random walk has constant mean.\n"""
    text = replace_once(text, old_mean, new_mean, "accelerometer-bias mean RW ablation")

    # 2) Residual accelerometer bias transition: remove OU contraction, retain
    #    the Level-2 coordinate rotation beta = R delta_b.
    old_transition = """        const Matrix3 Rw_step = ocean_imu::lie::Exp<T>(Vector3(X_.R * omega * Ts));\n        if constexpr (two_frame_bias) {\n            Phi.template block<3,3>(OFF_BG, OFF_BG) = Rw_step;\n        }\n        if constexpr (with_accel_bias) {\n            const T alpha_ba = acc_bias_updates_enabled_ ? std::exp(-Ts/tau_ba_) : T(1);\n            if constexpr (two_frame_bias)\n                Phi.template block<3,3>(OFF_BA, OFF_BA) = alpha_ba * Rw_step;\n            else\n                Phi.template block<3,3>(OFF_BA, OFF_BA) = alpha_ba * Matrix3::Identity();\n        }\n"""
    new_transition = """        const Matrix3 Rw_step = ocean_imu::lie::Exp<T>(Vector3(X_.R * omega * Ts));\n        if constexpr (two_frame_bias) {\n            Phi.template block<3,3>(OFF_BG, OFF_BG) = Rw_step;\n            if constexpr (with_accel_bias)\n                Phi.template block<3,3>(OFF_BA, OFF_BA) = Rw_step;\n        }\n        // Level 1 keeps Phi_ba,ba = I: historical random-walk model.\n"""
    text = replace_once(text, old_transition, new_transition,
                        "accelerometer-bias transition RW ablation")

    # 3) Gyro-bias Qd: restore the historical constant-rotation h^2/2,h^3/3
    #    approximation.  This is deliberately the code under test, not a claim
    #    that the approximation is exact.
    start = "        const Matrix3 R0 = X_.R;\n"
    end = "        Qd.template block<3,3>(OFF_BG, OFF_BG) = Bfinal * Q_bg_ * Bfinal.transpose() * Ts;\n"
    old_qd = """        const Matrix3 Rhat = X_.R;\n        const Matrix3 R_end =\n            two_frame_bias\n                ? Matrix3(Rhat * ocean_imu::lie::Exp<T>(\n                      Vector3((gyro_meas - X_.B.col(0)) * Ts)))\n                : Rhat;\n        const Matrix3 W  = Rhat * Q_gyro_ * Rhat.transpose();\n        const Matrix3 Wb = Rhat * Q_bg_   * Rhat.transpose();\n\n        const T h2 = Ts * Ts;\n        const T h3 = h2 * Ts;\n        const Matrix3 Q_bias_own = two_frame_bias\n            ? Matrix3(R_end * Q_bg_ * R_end.transpose()) : Q_bg_;\n        const Matrix3 Q_bias_cross = two_frame_bias\n            ? Matrix3(Rhat * Q_bg_ * R_end.transpose()) : Matrix3(Rhat * Q_bg_);\n\n        Qd.template block<3,3>(OFF_PHI, OFF_PHI) = W * Ts + Wb * (h3 / T(3));\n        Qd.template block<3,3>(OFF_PHI, OFF_BG) = -Q_bias_cross * (h2 / T(2));\n        Qd.template block<3,3>(OFF_BG, OFF_PHI) =\n            Qd.template block<3,3>(OFF_PHI, OFF_BG).transpose();\n        Qd.template block<3,3>(OFF_BG, OFF_BG) = Q_bias_own * Ts;\n"""
    text = replace_between(text, start, end, old_qd, "historical gyro-bias Qd ablation")

    # 4) Residual accelerometer bias process covariance: exact OU integral ->
    #    historical random-walk Q h.  Keep the hardened warm-up freeze gate so
    #    this arm changes the stochastic model, not startup policy.
    old_ba_q = """        if constexpr (with_accel_bias) {\n            if (acc_bias_updates_enabled_) {\n                const T qint = T(0.5) * tau_ba_ * (-std::expm1(T(-2)*Ts/tau_ba_));\n                Qd.template block<3,3>(OFF_BA, OFF_BA) =\n                    (two_frame_bias ? Matrix3(R_end * Q_ba_ * R_end.transpose()) : Q_ba_) * qint;\n            }\n        }\n"""
    new_ba_q = """        if constexpr (with_accel_bias) {\n            if (acc_bias_updates_enabled_) {\n                Qd.template block<3,3>(OFF_BA, OFF_BA) =\n                    (two_frame_bias ? Matrix3(R_end * Q_ba_ * R_end.transpose()) : Q_ba_) * Ts;\n            }\n        }\n"""
    text = replace_once(text, old_ba_q, new_ba_q,
                        "accelerometer-bias process RW ablation")
    return text


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
