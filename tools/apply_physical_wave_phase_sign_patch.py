#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_exact(path: Path, old: str, new: str, expected: int = 1) -> None:
    text = path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != expected:
        raise RuntimeError(f"{path}: expected {expected} occurrences, found {count}: {old!r}")
    path.write_text(text.replace(old, new), encoding="utf-8")


detector = ROOT / "src/wave_dir/WaveDirectionDetector.h"
replace_exact(
    detector,
    "    // default is correct for an up-positive vertical signal in which\n"
    "    // a_parallel and d(a_vertical)/dt have positive correlation for\n"
    "    // propagation along the positive axis.\n"
    "    Real convention_sign = Real(1);\n",
    "    // For linear deep-water orbital motion with up-positive vertical\n"
    "    // acceleration, a_parallel * d(a_vertical)/dt is negative for\n"
    "    // propagation along the positive axis.  Negating the coherence maps\n"
    "    // that physical convention to FORWARD.\n"
    "    Real convention_sign = Real(-1);\n")

main_test = ROOT / "tests/wave_dir/wave-direction-test.cpp"
replace_exact(main_test, "  cfg.convention_sign = 1.0f;\n",
              "  cfg.convention_sign = -1.0f;\n")
replace_exact(
    main_test,
    "    const float horizontal = float(true_sense) * 0.65f * std::cos(phase);\n"
    "    const float vertical_up = 0.50f * std::sin(phase);\n",
    "    // Linear deep-water orbital acceleration at a fixed observer:\n"
    "    // positive-axis propagation has negative correlation between\n"
    "    // a_parallel and d(a_up)/dt.\n"
    "    const float horizontal = -float(true_sense) * 0.65f * std::cos(phase);\n"
    "    const float vertical_up = 0.50f * std::sin(phase);\n")
replace_exact(
    main_test,
    "        const float horizontal = float(sense) * 0.70f * std::cos(phase);\n"
    "        const float vertical_up = 0.52f * std::sin(phase);\n",
    "        const float horizontal = -float(sense) * 0.70f * std::cos(phase);\n"
    "        const float vertical_up = 0.52f * std::sin(phase);\n")

iq_test = ROOT / "tests/wave_dir/wave-direction-iq-test.cpp"
replace_exact(
    iq_test,
    "  sense_config.min_axis_confidence = 20.0f;\n",
    "  sense_config.min_axis_confidence = 20.0f;\n"
    "  sense_config.convention_sign = -1.0f;\n")
replace_exact(
    iq_test,
    "    const float along = 0.65f * std::cos(phase);\n"
    "    const float vertical_up = 0.50f * std::sin(phase);\n",
    "    const float along = -0.65f * std::cos(phase);\n"
    "    const float vertical_up = 0.50f * std::sin(phase);\n")

print("Applied physical deep-water orbital phase sign convention.")
