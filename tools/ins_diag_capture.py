#!/usr/bin/env python3
"""Capture an AtomS3R INS diagnostic stream (SEA_STATE_SERIAL_DIAG=1) to a file.

    pip install pyserial
    python3 tools/ins_diag_capture.py --port /dev/ttyACM0 --out bench.log

Start the script, then reset or power-cycle the device so the capture begins at
boot: the USB serial port disappears during the reset and the script reopens
it, and the sketch waits up to five seconds for the reader before it starts.
Stop with Ctrl-C. Every received line is written unchanged; the script only
prints a status line every few seconds (frames per second, missing R lines,
dropped diagnostic lines reported by the device).

See docs/ins-device-diagnostics.md for what to record and how to replay it.
"""
from __future__ import annotations

import argparse
import sys
import time


def open_port(serial_module, port: str, baud: int):
    while True:
        try:
            return serial_module.Serial(port, baud, timeout=0.5)
        except (OSError, serial_module.SerialException):
            time.sleep(0.2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--port", required=True, help="serial device, e.g. /dev/ttyACM0 or COM5")
    parser.add_argument("--out", required=True, help="capture file to write (appended)")
    parser.add_argument("--baud", type=int, default=115200, help="ignored by USB CDC; default 115200")
    args = parser.parse_args()
    try:
        import serial  # pyserial
    except ImportError:
        print("pyserial is required: pip install pyserial", file=sys.stderr)
        return 2

    frames = gaps = 0
    expected = None
    last_report = time.monotonic()
    window_frames = 0
    dropped = "?"
    with open(args.out, "ab") as out:
        port = open_port(serial, args.port, args.baud)
        print(f"capturing {args.port} -> {args.out}; reset the device now to record from boot", file=sys.stderr)
        pending = b""
        try:
            while True:
                try:
                    chunk = port.read(4096)
                except (OSError, serial.SerialException):
                    port.close()
                    print("port lost; reopening", file=sys.stderr)
                    port = open_port(serial, args.port, args.baud)
                    pending = b""
                    continue
                if not chunk:
                    continue
                out.write(chunk)
                pending += chunk
                *lines, pending = pending.split(b"\n")
                for line in lines:
                    if line.startswith(b"R,"):
                        fields = line.split(b",")
                        try:
                            seq = int(fields[1], 16)
                        except (IndexError, ValueError):
                            continue
                        if seq == 0:
                            expected = 0
                            print("device boot (R sequence restarted)", file=sys.stderr)
                        if expected is not None and seq != expected:
                            gaps += 1
                        expected = seq + 1
                        frames += 1
                        window_frames += 1
                    elif line.startswith(b"D,"):
                        for field in line.split(b","):
                            if field.startswith(b"lines_dropped="):
                                dropped = field.split(b"=", 1)[1].decode(errors="replace")
                now = time.monotonic()
                if now - last_report >= 5.0:
                    rate = window_frames / (now - last_report)
                    print(f"frames {frames}  ({rate:.1f}/s)  R gaps {gaps}  device lines_dropped {dropped}",
                          file=sys.stderr)
                    out.flush()
                    last_report = now
                    window_frames = 0
        except KeyboardInterrupt:
            pass
    print(f"done: {frames} frames, {gaps} R gaps", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
