#!/usr/bin/env python3
"""Matched ordinary-aiding controls; finite engineering evidence, not a proof."""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
BACKENDS = (
    ("qMEKF", "ahrs", "rapid_motion-qmekf"),
    ("Mahony", "ahrs", "rapid_motion-mahony"),
    ("OU-II", "kalman_ou_ii", "rapid_motion-test"),
    ("OU-III", "kalman_ou_iii", "rapid_motion-test"),
    ("TFG", "kalman_tfg", "rapid_motion-test"),
    ("NLO", "nlo", "rapid_motion-test"),
    ("PII", "pii_observer", "rapid_motion-test"),
)


def sha(path):
    """Return the SHA-256 hex digest of the bytes at path."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    """Build all backends, run matched controls, and write metrics with a source-bound manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eigen", default=os.environ.get("EIGEN_INCLUDE_DIR"))
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports/results/rapid_motion")
    args = parser.parse_args()
    rows = []
    recovery_rows = []
    binaries = {}
    regressions = {}
    for number, (name, directory, target) in enumerate(BACKENDS, 1):
        folder = ROOT / "tests" / directory
        # Bind the recorded source to a fresh executable, including headers
        # in the older test Makefiles that do not emit automatic dependencies.
        command = ["make", "-B", "-C", str(folder), target]
        if args.eigen:
            command.append("EIGEN_DIR=" + args.eigen)
        subprocess.run(command, check=True)
        binary = folder / target
        recovery = subprocess.run([str(binary), "--diagnose-recovery"], check=True,
                                  capture_output=True, text=True)
        cases = list(csv.DictReader(io.StringIO(recovery.stdout)))
        if len(cases) != 4 or any(int(r["backend"]) != number for r in cases):
            raise ValueError(f"{name}: incomplete inherited-state controls")
        for row in cases:
            row["backend"] = name
        for start in (0, 2):
            if any(cases[start][key] != cases[start+1][key] for key in
                   ("entry_attitude_deg", "entry_bias_rad_s")):
                raise ValueError(f"{name}: recovery pair inherited different states")
        recovery_rows.extend(cases)
        for mode in ("--diagnose", "--diagnose-force-equation"):
            replay = subprocess.run([str(binary), mode], check=True,
                                    capture_output=True, text=True)
            cases = list(csv.DictReader(io.StringIO(replay.stdout)))
            if len(cases) != 24 or any(int(r["backend"]) != number or
                                       int(r["applied_recoveries"]) != 0 for r in cases):
                raise ValueError(f"{name}: incomplete or protected causal controls")
            for row in cases:
                row["backend"] = name
                row["force_equation_control"] = int(mode.endswith("force-equation"))
            rows.extend(cases)
        regression = subprocess.run([str(binary)], check=True,
                                    capture_output=True, text=True)
        records = list(csv.reader(io.StringIO(regression.stdout)))
        rapid = [r for r in records if r[0] == "rapid" and r[2] == "1"]
        capture = [r for r in records if r[0] == "capture"]
        antipodes = [r for r in records if r[0] == "antipode"]
        dropouts = [r for r in records if r[0] == "accel-dropout"]
        if (len(rapid), len(capture), len(antipodes), len(dropouts)) != (22, 14, 2, 3):
            raise ValueError(f"{name}: incomplete protection regressions")
        regressions[name] = {
            "protected_rapid_cases": len(rapid),
            "capture_cases": len(capture),
            "dropout_cases": len(dropouts),
            "maximum_rapid_attitude_deg": max(float(r[8]) for r in rapid),
            "maximum_rapid_tail_deg": max(float(r[9]) for r in rapid),
            "maximum_capture_tail_deg": max(float(r[4]) for r in capture),
            "local_antipode_tail_deg": float(antipodes[0][3]),
            "protected_antipode_tail_deg": float(antipodes[1][3]),
        }
        binaries[str(binary.relative_to(ROOT))] = sha(binary)
    args.output.mkdir(parents=True, exist_ok=True)
    output = args.output / "causal-controls.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    sources = sorted((ROOT / "src").rglob("*.h")) + [
        ROOT / "tests/common/RapidMotionReplay.h",
        ROOT / "tests/common/RapidMotionBackends.h",
        Path(__file__).resolve(),
    ]
    recovery_output = args.output / "recovery-controls.csv"
    with recovery_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(recovery_rows[0]))
        writer.writeheader()
        writer.writerows(recovery_rows)
    manifest = {
        "inherited_state_recovery": {
            "scope": "Ordinary aiding through the 720 deg/s, 15 cm turn; paired "
                     "continuations retain every state. Protection is enabled "
                     "only after t=11 s in the late-protection arm. No reseeding. "
                     "This deliberately bypasses prevention during the turn.",
            "cases": len(recovery_rows),
            "csv_sha256": sha(recovery_output),
            "tail_seconds_after_turn": 7,
            "diagnostic_target_deg": 5,
            "late_protection_cases_above_target": sum(
                int(r["protection_after_turn"]) == 1 and
                float(r["tail_attitude_deg"]) >= 5 for r in recovery_rows),
            "qualification": "Diagnostic, not a recovery acceptance gate.",
        },
        "scope": "Production cores and matched delivery adapters; ordinary aiding "
                 "enabled, no protection/reacquisition applied. Synthetic motion, "
                 "not complete device execution or hardware qualification.",
        "cases_per_backend": 48,
        "force_equation_control": "Subtract exact alpha cross r + omega cross "
                                  "(omega cross r) from the accelerometer input. "
                                  "Known-truth causal control, not a shipping "
                                  "algorithm for an unknown moving pivot.",
        "tail_seconds_after_turn": 7,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
        "binary_sha256": binaries,
        "protection_regressions": regressions,
        "csv_sha256": sha(output),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(rows)} matched cases to {output}")


if __name__ == "__main__":
    main()
