#!/usr/bin/env python3
"""Copyright 2026, Mikhail Grushinskiy. Compile the device wizard's config references."""
import os
from pathlib import Path
import re
import shlex
import subprocess

root = Path(__file__).resolve().parents[2]
wizard = root / "src/AtomS3R/AtomS3R_ImuCalWizard.h"
text = wizard.read_text()
definition = re.search(r"^struct ImuCalWizardCfg \{.*?^\};", text, re.M | re.S)
assert definition, "Wizard configuration declaration not found"
members = sorted(set(re.findall(r"ImuCalWizardCfg::(\w+)", text)))
assert members, "No wizard configuration uses found"

# Host tests do not include the Arduino/FreeRTOS UI. Compile the actual config
# declaration and every member used by that UI, with only its platform type
# supplied here. This checks name resolution, not the full MCU build.
source = "#include <cstdint>\nusing StackType_t = uint32_t;\n"
source += definition.group() + "\nvoid check_wizard_config() {\n"
source += "".join(f"  (void)ImuCalWizardCfg::{member};\n" for member in members)
source += "}\n"
subprocess.run(
    shlex.split(os.environ.get("CXX", "g++"))
    + ["-std=c++17", "-x", "c++", "-fsyntax-only", "-"],
    input=source,
    text=True,
    check=True,
)
print(f"test_wizard_config: {len(members)} device configuration members compile")
