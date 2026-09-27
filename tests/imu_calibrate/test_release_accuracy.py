#!/usr/bin/env python3
"""Keep the deployed accelerometer procedure at least as accurate as v2.3.2.

Reference: 2ba31e4e255a2b03f014cfbb9c9352d58d49f791, its own 30-seed
accel_cal-test campaign, RICH+NEW (the procedure actually deployed in that tag).
The OLD column of that campaign predates v2.3.2 and is not this baseline.
"""
import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent
with (root / 'calibrate_v232_accel_reference.csv').open() as f:
    reference = list(csv.DictReader(f))
with (root / 'calibrate_accel_campaign_summary.csv').open() as f:
    current = {r['scenario']: r for r in csv.DictReader(f) if r['method'] == 'RICH+NEW'}
for old in reference:
    new = current[old['scenario']]
    assert float(new['success_rate']) >= float(old['success_rate']), old['scenario']
    for key in old.keys() - {'scenario', 'success_rate'}:
        a, b = float(old[key]), float(new[key])
        # Allow report rounding and small compiler arithmetic differences.
        assert math.isfinite(b) and b <= a * 1.02 + 1e-5, (old['scenario'], key, a, b)
print(f'v2.3.2 accelerometer reference: {len(reference)} scenarios PASS')
