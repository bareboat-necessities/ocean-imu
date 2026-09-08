"""Insert passive runtime callbacks; stripping them recovers shipping bytes.

The ordinary simulator and harness are reused, including their noise, startup,
25 Hz magnetic clock and quality gates. No shipping file is edited. A dual run
must additionally recover all ordinary simulator output byte for byte.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

SIM = "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"
MEKF = "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
FUSION = "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"


def instrument(original, path):
    source, edits = original, []

    def add(anchor, call, before=True):
        nonlocal source
        if source.count(anchor) != 1:
            raise ValueError(f"nonunique runtime trace anchor in {path}: {anchor}")
        statement = "\n" + call + " // OU3_BRMM_READ_ONLY\n"
        source = source.replace(anchor, statement+anchor if before else anchor+statement)
        edits.append(statement)

    if path == SIM:
        add('#include "util/W3dSimCommon.h"', '#include "ou3_brmm_runtime_trace.h"', False)
        add("        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);",
            "        ou3_brmm_trace::begin(fusion_, dt);")
        add("        fusion_.updateMag(mag_body_ned);", "        ++ou3_brmm_trace::get().mag_delivered;")
        # snapshot follows both IMU and asynchronous magnetic processing in the
        # unchanged harness. Logging at update() return would misalign mag events.
        add("    FilterSnapshot snapshot() const override {", "        ou3_brmm_trace::finish(fusion_);", False)
    elif path == MEKF:
        add("    last_dt_ = Ts;   // Remember last dt", "    ++ou3_brmm_trace::get().predict;")
        add("    last_acc_diag_ = MeasDiag3{};", "    ++ou3_brmm_trace::get().acc_calls;")
        add("    last_mag_diag_ = MeasDiag3{};", "    ++ou3_brmm_trace::get().mag_calls;")
        add("    xext.noalias() += K * r;          // State update",
            '    ou3_brmm_trace::measurement("acc", f_cog_b, R_wb(), Racc);')
        add("    // State + covariance update\n    xext.noalias() += K * r;",
            '    ou3_brmm_trace::measurement("mag", v2hat, R_wb(), Rmag);')
        add("    constexpr int off_S = OFF_S;   // offset of S block (3 states)",
            "    ++ou3_brmm_trace::get().s_calls;")
        add("    xext.noalias() += K * r;            // State update",
            '    ou3_brmm_trace::measurement("S", r, Matrix3::Identity(), R_S);')
        add("    const Vector3 dtheta = xext.template segment<3>(0);",
            "    ++ou3_brmm_trace::get().resets;")
        add("    void set_mag_world_ref(const Vector3& B_world) {",
            "        ++ou3_brmm_trace::get().reference_writes;", False)
    elif path == FUSION:
        add("        const Eigen::Vector3f acc_in = accel_guard_.step(acc, dt);",
            "        ou3_brmm_trace::guard(acc, acc_in, accel_guard_.engagement(), accel_guard_.excessRms());", False)
        add("            if (tilt_over_limit_sec_ >= TILT_RESET_HOLD_SEC && tilt_reset_cooldown_sec_ <= 0.0f) {",
            "                ++ou3_brmm_trace::get().tilt_resets;", False)
    else:
        raise ValueError("unknown runtime overlay path")
    stripped = source
    for edit in edits:
        if stripped.count(edit) != 1:
            raise ValueError("runtime trace insertion identity lost")
        stripped = stripped.replace(edit, "")
    if stripped != original:
        raise ValueError("runtime overlay changed shipping arithmetic")
    return source, edits


def build(repository, output):
    manifest = {"stripping_recovers_shipping_bytes": True,
                "runtime_transparency_requires_dual_run": True, "files": {}}
    for path in (SIM, MEKF, FUSION):
        original = (repository/path).read_text()
        generated, edits = instrument(original, path)
        target = output/path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(generated)
        manifest["files"][path] = {
            "shipping_sha256": hashlib.sha256(original.encode()).hexdigest(),
            "overlay_sha256": hashlib.sha256(generated.encode()).hexdigest(),
            "insertions": edits}
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.repository, args.output)
