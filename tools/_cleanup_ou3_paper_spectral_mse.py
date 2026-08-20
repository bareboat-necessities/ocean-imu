#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "doc" / "kalman_ou_iii"


def patch(name, replacements):
    path = DOC / name
    text = path.read_text(encoding="utf-8")
    for old, new in replacements:
        if old not in text:
            raise SystemExit(f"{name}: missing expected text: {old!r}")
        text = text.replace(old, new)
    path.write_text(text, encoding="utf-8")


patch("w3d-reduced-mse-validation.tex-part", [
    ("direct paired comparison against the applied law", "direct paired comparison against the legacy amplitude-coupled cubic comparator"),
    ("measures the resulting law against the applied schedule", "measures the resulting law against the legacy cubic comparator"),
    ("at the applied operating points", "at the calibrated operating points"),
    ("Writing the applied\nform as", "Writing the deployed SpectralMSE\nform as"),
    ("in the applied parameterization", "in the legacy cubic parameterization"),
    ("both\nthe applied $p_a=1$", "both\nthe legacy cubic $p_a=1$"),
    ("the\napplied law gives $3.61$", "the\nlegacy cubic law gives $3.61$"),
    ("compared against the applied law under", "compared against the legacy cubic comparator under"),
    ("against the\napplied schedule. Vertical-displacement", "against the\nlegacy cubic comparator. Vertical-displacement"),
    ("fitted to the\napplied schedule's own worst values", "fitted to the\nlegacy cubic schedule's own worst values"),
])

patch("w3d-reduced-mse-envelope.tex-part", [
    ("Neither exponent is fitted to the applied adaptation law.",
     "Neither exponent was fitted to the legacy cubic comparator; they are the exponents of the deployed SpectralMSE schedule."),
    ("MSE/applied $-1$", "MSE/cubic $-1$"),
])

patch("w3d-adaptation-coefficient-investigation.tex-part", [
    ("followed by the same cadence normalization as the applied filter.",
     "followed by the same reference-cadence normalization as the legacy cubic comparator."),
])

# Assert that the primary implementation narrative no longer presents the
# former cubic schedule as current.
checks = {
    "w3d-adaptation-motivation.tex-part": [
        "The applied filter retains the dimensionally cubic base",
        "The resulting applied schedule is summarized by",
    ],
    "w3d-fus-methods.tex-part": [
        "r_{S,\\star}^{\\mathrm{base}}=c_R\\sigma_{aw,\\star}\\tau_\\star^3",
        "adaptation commit / EMA horizon & $0.1/1.8$ s",
    ],
    "w3d-reduced-mse-validation.tex-part": [
        "against the applied law",
        "against the applied schedule",
        "in the applied parameterization",
        "the applied law gives $3.61$",
    ],
    "w3d-ou-robustness.tex-part": [
        "according to the applied relation",
    ],
}
for name, forbidden in checks.items():
    text = (DOC / name).read_text(encoding="utf-8")
    for phrase in forbidden:
        if phrase in text:
            raise SystemExit(f"{name}: stale deployed-cubic wording remains: {phrase}")

print("OU-III manuscript terminology cleanup complete")
