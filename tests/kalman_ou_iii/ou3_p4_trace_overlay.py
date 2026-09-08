"""Generate a read-only event-trace overlay of the unchanged shipping header.

Only the enumerated callback statements are inserted. Removing them must
recover the shipping header byte-for-byte. CI additionally compares the
instrumented observer's complete existing capture with the uninstrumented
binary. This is instrumentation, not a second transition implementation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def instrument(source):
    edits = []

    def insert(anchor, statement, *, before=True):
        nonlocal source
        if source.count(anchor) != 1:
            raise ValueError("shipping trace anchor is not unique: " + anchor[:100])
        addition = "\n    " + statement + " // OU3_READ_ONLY_TRACE\n"
        source = source.replace(anchor, addition+anchor if before else anchor+addition)
        edits.append(addition)

    insert("    last_dt_ = Ts;   // Remember last dt",
           'ou3_source_trace::state("prediction_enter", *this);')
    insert("    last_gyr_bias_corrected = gyr - gyro_bias;",
           'ou3_source_trace::gyro(gyr);', before=False)
    insert("    apply_pending_aw_covariance_inflation_();\n    symmetrize_Pext_();   // Symmetry hygiene",
           'ou3_source_trace::state("prediction", *this);')
    insert("    apply_pending_aw_covariance_inflation_();\n    symmetrize_Pext_();   // Symmetry hygiene",
           'ou3_source_trace::state("aw_floor", *this);', before=False)
    insert("    xext.noalias() += K * r;          // State update",
           'ou3_source_trace::measurement("accelerometer", *this, r, acc_meas, Racc, k_a_ * (tempC - tempC_ref));')
    insert("    // State + covariance update\n    xext.noalias() += K * r;",
           'ou3_source_trace::measurement("magnetometer", *this, r, mag_meas, Rmag, Vector3::Zero());')
    insert("    xext.noalias() += K * r;            // State update",
           'ou3_source_trace::measurement("S_zero", *this, r, Vector3::Zero(), R_S, Vector3::Zero());')
    insert("    const Vector3 dtheta = xext.template segment<3>(0);",
           'ou3_source_trace::state("injected", *this);')
    insert("    project_acc_bias_();",
           'ou3_source_trace::state("reset", *this);')
    insert("    project_acc_bias_();",
           'ou3_source_trace::state("projection", *this);', before=False)
    return source, edits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    relative = Path("kalman_ou_iii/Kalman3D_Wave_OU_III.h")
    path = args.repository / "src" / relative
    original = path.read_text()
    generated, edits = instrument(original)
    stripped = generated
    for edit in edits:
        if stripped.count(edit) != 1:
            raise ValueError("trace statement identity lost")
        stripped = stripped.replace(edit, "")
    if stripped != original:
        raise ValueError("overlay changed shipping arithmetic")
    target = args.output / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(generated)
    manifest = {"shipping_header_sha256": hashlib.sha256(original.encode()).hexdigest(),
                "overlay_header_sha256": hashlib.sha256(generated.encode()).hexdigest(),
                "read_only_insertions": edits, "stripping_recovers_shipping_bytes": True,
                "runtime_bitwise_transparency_still_requires_dual_run": True}
    (args.output / "trace-manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    main()
