"""Close the full-state process-envelope step used by the analytical proof."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 60
D = Decimal
ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou_live_basin_interval_proof.py"


class ComputerAssistedProcessBlockTests(unittest.TestCase):
    def test_reported_translational_window_encloses_other_default_process_blocks(self):
        cp = subprocess.run(
            [sys.executable, str(TOOL), "--repo-root", str(ROOT), "--json"],
            check=True,
            text=True,
            capture_output=True,
        )
        v = {k: D(x) for k, x in json.loads(cp.stdout).items()}
        T = v["horizon_max_s"]

        # Scaled attitude/gyro-bias process envelope. Rotation preserves the
        # gyro-noise norm. Bias driving noise contributes directly to b_g and
        # through the attitude integral B(t), with ||B(t)|| <= t.
        qg = (D("0.00157") / D("0.087")) ** 2
        qbg = D("1e-11") / D("0.001") ** 2
        cb = D("0.001") / D("0.087")
        w_att = D(3) * (
            qg * T + qbg * (T + cb * cb * T**3 / D(3))
        )

        # Scaled residual-bias covariance never exceeds its OU stationary
        # covariance in the comparison estimator used by the proof.
        pba = max(
            (D("0.004") / D("0.5")) ** 2,
            D("2.5e-7") * D("5000") / (D(2) * D("0.5") ** 2),
        )
        w_ba = D(3) * pba

        self.assertGreater(v["process_window_upper"], w_att)
        self.assertGreater(v["process_window_upper"], w_ba)

        # The finite-horizon translational controllability constant is also a
        # valid full-state process floor because the other two process blocks
        # have much larger one-step floors.  For attitude/gyro bias, use
        # ||B(t)|| <= cb*t in scaled coordinates and spend at most half of the
        # direct gyro-noise term on the cross term.
        h = v["h_min"]
        hmax = v["h_max"]
        eta = min(
            D(1),
            D(3) * qg / (D(2) * qbg * cb * cb * hmax * hmax),
        )
        qatt_step = h * min(qg / D(2), qbg * eta / (D(1) + eta))

        # The residual accelerometer-bias OU block has exact O(h) diffusion.
        qba_step = v["qba_step_lower"]

        self.assertGreater(qatt_step, v["alpha_controllability_lower"])
        self.assertGreater(qba_step, v["alpha_controllability_lower"])
        self.assertGreater(qatt_step, D("1e-9"))
        self.assertGreater(qba_step, D("1e-9"))

        # The old final-sample translational h^7 bound remains much smaller and
        # is retained only as a conservatism diagnostic.
        self.assertGreater(qatt_step, v["q_step_lower"])
        self.assertGreater(qba_step, v["q_step_lower"])


if __name__ == "__main__":
    unittest.main()
