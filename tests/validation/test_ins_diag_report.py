"""tools/ins_diag_report.py finds the device events it is meant to report."""
from __future__ import annotations

import io
import struct
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ins_diag_report as report  # noqa: E402


def bits(x: float) -> str:
    return format(struct.unpack("<I", struct.pack("<f", x))[0], "x")


def r_line(seq: int, us: int, gyro: float = 0.0, mag_seq: int = 1, mag_us: int = 0) -> str:
    floats = [0.0, 0.0, -9.8, gyro, 0.0, 0.0, 20.0, 0.0, 45.0, 30.0]
    timing = [1, us // 1000, mag_seq, mag_us, mag_us, mag_us + 300]
    return ",".join(["R", format(seq, "x"), format(us, "x"), "3"] + [bits(f) for f in floats] +
                    [format(v, "x") for v in timing])


def d_line(t_ms: int, **kv) -> str:
    base = {"t_ms": t_ms, "frames": 10, "overruns": 0, "tick_max_us": 900, "dt_max_us": 10000,
            "lines_dropped": 0, "live": 1, "north": 1, "reseeds": 0, "roll": 1.0, "pitch": 0.0,
            "hdg": 90.0, "hdg_cmp": 91.0, "hdg_ok": 1, "heave": 0.0, "env": 0.0,
            "gate_trusted": 1, "gyro_fast": 0, "gate_good_s": 2.0, "mag_tries": 2, "mag_applied": 2}
    base.update(kv)
    return "D," + ",".join(f"{k}={v}" for k, v in base.items())


class InsDiagReportTest(unittest.TestCase):
    def build(self) -> list[str]:
        lines = ["E,reset_fusion,seq=0"]
        us = 2_000_000
        for seq in range(0, 3000):
            if seq == 1500:
                continue  # one lost R line
            gap = 60_000 if seq == 2000 else 10_000  # one long frame gap
            us += gap
            lines.append(r_line(seq, us, gyro=0.8 if seq % 2 else 0.0, mag_seq=1 + seq // 4,
                                mag_us=us - 1000))
            if seq == 0:
                lines += ["H,sketch=atomS3R_ins_kalman_ou3,format=1,loop_hz=200.0,seq=1",
                          "I,bmi270=1,acc_conf=0xa8,acc_odr_hz=100.00,gyr_conf=0xa9,gyr_odr_hz=200.00,"
                          "aux_odr_hz=25.00,bmm150=1,bmm_odr_hz=30.0,bmm_rep_xy=23,bmm_rep_z=40,mag_source=1",
                          "C,1,1," + "00" * 372]
            if seq % 10 == 9:
                t_ms = (us - 0) // 1000
                k = seq // 10
                kv = {}
                if k < 20:
                    kv["live"] = 0
                if 100 <= k < 130:
                    kv["hdg"] = 270.0                    # north/south-scale disagreement
                if 150 <= k < 200:
                    kv["heave"] = 0.25                   # frozen heave
                else:
                    kv["heave"] = 0.001 * (k % 7)
                if k >= 220:
                    kv["reseeds"] = 1
                if 230 <= k < 240:
                    kv["roll"] = 75.0
                lines.append(d_line(t_ms, **kv))
        return lines

    def test_report_lists_each_event(self):
        boots = report.parse(self.build())
        self.assertEqual(len(boots), 1)
        out = io.StringIO()
        report.report_boot(0, boots[0], out)
        text = out.getvalue()
        self.assertIn("acc_odr 100.00 Hz", text)
        self.assertIn("missing R lines 1", text)
        self.assertIn(">50 ms: 1", text)
        self.assertIn("frames above 30 dps 50.0%", text)
        self.assertIn("live: at", text)
        self.assertIn("NORTH/SOUTH-SCALE disagreement", text)
        self.assertIn("heave unchanged for >= 3 s", text)
        self.assertIn("re-seeds (yaw sigma jump > 30 deg): 1", text)
        self.assertIn("large tilt", text)
        self.assertIn("reset_fusion", text)

    def test_fifo_timing_lines(self):
        lines = [r_line(0, 1000), r_line(1, 6000),
                 "F,ok=1,failure=none,acc_conf=0xa9,gyr_conf=0xa9,aux_conf=0x06,aux_offset_ms=40.0,mag_cal=1",
                 "X,ok=1,delay_ms=28.0,width_ms=2.0,reason=Delay validated",
                 "X,ok=0,delay_ms=0.0,width_ms=0.0,reason=Weak timing signal",
                 "P,window_s=120.0,rate_hz=199.600,lost=0,stalls=12,stall_max_ms=1500,fifo_max_bytes=4200,"
                 "spacing_min_us=4995.0,spacing_max_us=5015.0,clock_ppm=2000,aux_hz=25.000,drdy_fail=0,"
                 "mag_repeat=3,aux_interval_bad=0,delay_fits=4,delay_ms=28.5,delay_spread_ms=2.0,"
                 "rate_ok=1,no_loss_ok=1,spacing_ok=1,mag_ok=1,delay_ok=1,all_ok=1"]
        boots = report.parse(lines)
        out = io.StringIO()
        report.report_boot(0, boots[0], out)
        text = out.getvalue()
        self.assertIn("FIFO acquisition setup: ok", text)
        self.assertIn("overall  PASS", text)
        self.assertIn("delay fits: 1 of 2 valid, 28 ms; rejected: Weak timing signal", text)

    def test_replay_output_without_frames(self):
        boots = report.parse([d_line(1000), d_line(1100)])
        self.assertEqual(len(boots), 1)
        self.assertEqual(len(boots[0].state), 2)

    def test_new_boot_starts_new_segment(self):
        lines = [r_line(0, 1000), r_line(1, 11000), "E,reset_fusion,seq=0", r_line(0, 5000), r_line(1, 15000)]
        boots = report.parse(lines)
        self.assertEqual([len(b.seqs) for b in boots], [2, 2])
        self.assertEqual(boots[1].events, ["reset_fusion,seq=0"])


if __name__ == "__main__":
    unittest.main()
