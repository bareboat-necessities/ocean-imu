#!/usr/bin/env python3
"""One-shot branch migration for independent TFG X/Y R_S coefficients.

This helper is used only to materialize the source edit on the experiment branch;
it is removed again after the generated source commit exists.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{label}: expected one match, found {count}")
    return text.replace(old, new, 1)


header = ROOT / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"
text = header.read_text()
text = replace_once(
    text,
    "    void setRSXYFactor(float k)    { if (k > 0.0f && std::isfinite(k)) R_S_xy_factor_ = k; }",
    "    void setRSXFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_x_factor_ = k; }\n"
    "    void setRSYFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_y_factor_ = k; }\n"
    "    [[nodiscard]] float getRSXFactor() const noexcept { return R_S_x_factor_; }\n"
    "    [[nodiscard]] float getRSYFactor() const noexcept { return R_S_y_factor_; }",
    "setter API",
)
text = replace_once(
    text,
    "        mekf_.set_RS_noise(Vector3f(rs * R_S_xy_factor_, rs * R_S_xy_factor_, rs));",
    "        mekf_.set_RS_noise(Vector3f(rs * R_S_x_factor_, rs * R_S_y_factor_, rs));",
    "R_S application",
)
text = replace_once(
    text,
    "    // TFG-specific physical prior/anisotropy coefficients remain the values\n"
    "    // independently measured for TFG.  Front-end statistical coefficients are\n"
    "    // what are kept at OU-III parity.\n"
    "    float tau_coeff_ = 1.0f;\n"
    "    float sigma_coeff_ = 1.0f;\n"
    "    float R_S_coeff_ = 0.28f;\n"
    "    float S_factor_ = 1.00f;\n"
    "    float R_S_xy_factor_ = 1.15f;",
    "    // TFG-specific physical OU prior coefficients remain independently fitted.\n"
    "    // Integral-state regularization is isotropic by default, matching OU-III;\n"
    "    // X/Y factors are independent opt-in experiment/tuning knobs.\n"
    "    float tau_coeff_ = 1.0f;\n"
    "    float sigma_coeff_ = 1.0f;\n"
    "    float R_S_coeff_ = 0.28f;\n"
    "    float S_factor_ = 1.00f;\n"
    "    float R_S_x_factor_ = 1.0f;\n"
    "    float R_S_y_factor_ = 1.0f;",
    "default factors",
)
header.write_text(text)

sim = ROOT / "tests/kalman_tfg/kalman_tfg-sim.cpp"
text = sim.read_text()
text = replace_once(
    text,
    "        if (env_float(\"TFG_R_S_XY_FACTOR\", v))  fusion_.setRSXYFactor(v);",
    "        if (env_float(\"TFG_R_S_X_FACTOR\", v))   fusion_.setRSXFactor(v);\n"
    "        if (env_float(\"TFG_R_S_Y_FACTOR\", v))   fusion_.setRSYFactor(v);",
    "simulator overrides",
)
sim.write_text(text)

print("materialized independent TFG R_S X/Y factors; defaults are 1.0/1.0")
