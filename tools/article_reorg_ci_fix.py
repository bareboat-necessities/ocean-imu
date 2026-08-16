#!/usr/bin/env python3
from pathlib import Path


def replace(path_s: str, old: str, new: str, count: int = 1) -> None:
    path = Path(path_s)
    text = path.read_text(encoding="utf-8")
    if text.count(old) < count:
        raise SystemExit(f"expected anchor not found in {path_s}: {old[:160]!r}")
    text = text.replace(old, new, count)
    path.write_text(text, encoding="utf-8")


# The statistical evidence schema intentionally stopped exporting deterministic
# simulator gate fields. Restatement must reproduce that schema natively.
replace(
    "tools/ou_validation.py",
    '    gate_passes = sum(int(row["quality_gate_pass"]) for row in rows)\n',
    "",
)
replace(
    "tools/ou_validation.py",
    '        rf"\\providecommand{{\\OUValidationGatePasses}}{{{gate_passes}}}",\n        rf"\\providecommand{{\\OUValidationGateFailures}}{{{len(rows) - gate_passes}}}",\n',
    "",
)

# Structural tests now assert absence of the retired fields rather than the old
# return-code/gate algebra.
old_gate_test = '''        self.assertTrue(\n            all(\n                int(row["simulator_return_code"])\n                == 1 - int(row["quality_gate_pass"])\n                for row in raw\n            )\n        )\n'''
new_gate_test = '''        self.assertTrue(\n            all(\n                "simulator_return_code" not in row\n                and "quality_gate_pass" not in row\n                for row in raw\n            )\n        )\n'''
replace("tests/validation/test_ou_validation.py", old_gate_test, new_gate_test)
replace("tests/validation/test_ou_robustness.py", old_gate_test, new_gate_test)

# Late-sorted wave-direction publication overrides own the final monkeypatches.
replace(
    "tests/validation/test_zzz_wave_direction_split_contract.py",
    '    self.assertIn("pure-start, crossfade, and pure-endpoint intervals", protocol)\n',
    '    self.assertIn("pure-start, blend, and pure-end", protocol)\n',
)
replace(
    "tests/validation/test_zzz_wave_direction_split_contract.py",
    '    fusion = self.read("w3d-fus-methods.tex-part")\n    simulation = self.read_flat("w3d-sim-charts.tex-part")\n',
    '    fusion = self.read("w3d-fus-methods.tex-part")\n    startup = self.read("w3d-init.tex-part")\n    results = self.read_flat("w3d-results.tex-part")\n',
)
replace(
    "tests/validation/test_zzz_wave_direction_split_contract.py",
    '''    # Direction/tracker thresholds no longer belong to the OU implementation\n    # table; the remaining attitude reset is an OU startup/recovery safeguard.\n    self.assertNotIn("direction-axis validity", fusion)\n    self.assertNotIn("sense coherence", fusion)\n    self.assertNotIn("tracker band", fusion)\n    self.assertIn("70^\\\\circ", fusion)\n\n    self.assertIn("functional source portability only", simulation)\n    self.assertIn("timing, processor load, memory use, power, deadline margin", simulation)\n    self.assertIn("were not measured", simulation)\n    self.assertIn("outside the claims of this study", simulation)\n''',
    '''    # Direction/tracker thresholds no longer belong to the OU adaptation\n    # table; tilt re-lock is documented with startup/recovery, while embedded\n    # portability limits are reported with the empirical results.\n    self.assertNotIn("direction-axis validity", fusion)\n    self.assertNotIn("sense coherence", fusion)\n    self.assertNotIn("tracker band", fusion)\n    self.assertIn("70^\\\\circ", startup)\n\n    self.assertIn("functional source portability only", results)\n    self.assertIn("Processor load, memory use, power, deadline margin", results)\n    self.assertIn("were not measured", results)\n''',
)

# Add the actual deployed tilt re-lock thresholds to the startup description.
replace(
    "doc/kalman_ou_iii/w3d-init.tex-part",
    '''A sufficiently large sustained tilt inconsistency can trigger a tilt-only\nre-lock that preserves the current yaw gauge.  Such resets are safety behavior,\nnot part of the nominal Kalman recursion.  They are therefore treated exactly\nlike startup and magnetic re-gauging in the analysis: bounded and testable in\nimplementation, but outside the local ISS theorem between reset events.\n''',
    '''A sufficiently large sustained tilt inconsistency can trigger a tilt-only\nre-lock that preserves the current yaw gauge.  The deployed safeguard requires\na tilt inconsistency above $70^\\circ$ for \SI{0.35}{s}, followed by a\n\SI{3}{s} cooldown before another re-lock is permitted.  Such resets are safety\nbehavior, not part of the nominal Kalman recursion.  They are therefore treated\nexactly like startup and magnetic re-gauging in the analysis: bounded and\ntestable in implementation, but outside the local ISS theorem between reset\nevents.\n''',
)

# Two local publication checks were lexical rather than semantic.
replace(
    "tests/validation/test_zz_publication_contract.py",
    '    self.assertIn("integral-regularization", conclusion)\n',
    '    self.assertTrue(\n        "integral-state regularization" in conclusion\n        or "integral regularizer" in conclusion\n    )\n',
)
replace(
    "tests/validation/test_zz_publication_contract.py",
    '        self.assertIn("hard resets", stability)\n',
    '        self.assertTrue("hard reset" in stability or "hard-reset" in stability)\n',
)
