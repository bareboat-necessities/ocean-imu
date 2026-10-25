#!/bin/bash
# ESP32-S3 size/stack probe of the three-axis CARMA filters (see
# carma3_size_probe.cpp).  Usage: build.sh <xtensa toolchain bin dir> [out dir]
# The toolchain is Espressif's crosstool-NG release (esp-14.2.0_20241119, the
# GCC of arduino-esp32 3.3.x); flags follow arduino-esp32's C++ defaults.
set -euo pipefail
TC=${1:?xtensa toolchain bin dir}
OUT=${2:-$(mktemp -d)}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../.." && pwd)
EIGEN=${EIGEN_DIR:-/usr/include/eigen3}
for v in 0 1 2 3 4 5 6; do
  (cd "$OUT" && "$TC/xtensa-esp32s3-elf-g++" -std=gnu++2a -Os -mlongcalls -ffunction-sections -fdata-sections \
      -fexceptions -fno-rtti -fstack-usage -DCARMA3_PROBE_VARIANT=$v -I"$ROOT/src" -isystem "$EIGEN" \
      -c "$HERE/carma3_size_probe.cpp" -o "$OUT/v$v.o")
  "$TC/xtensa-esp32s3-elf-objdump" -h "$OUT/v$v.o" > "$OUT/v$v.h.txt"
  for s in fusion mekf; do
    "$TC/xtensa-esp32s3-elf-objdump" -s -j ".rodata.carma3_probe_sizeof_$s" "$OUT/v$v.o" > "$OUT/v$v.$s.txt"
  done
done
python3 - "$OUT" <<'PY'
import csv, re, sys
from pathlib import Path
out = Path(sys.argv[1])
names = ["OU3", "Z41", "Z52", "A41", "A52", "M45", "A52S"]
wr = csv.writer(sys.stdout, lineterminator="\n")
wr.writerow(["variant", "code_bytes", "rodata_bytes", "data_bytes", "bss_bytes", "sizeof_fusion", "sizeof_filter",
             "max_stack_frame", "max_stack_function"])
for v, n in enumerate(names):
    tot = {"code": 0, "rodata": 0, "data": 0, "bss": 0}
    for line in (out / f"v{v}.h.txt").read_text().splitlines():
        f = line.split()
        if len(f) < 3 or not f[0].isdigit():
            continue
        sec, size = f[1], int(f[2], 16)
        if sec.startswith((".text", ".literal")): tot["code"] += size
        elif sec.startswith(".rodata"): tot["rodata"] += size
        elif sec.startswith(".data"): tot["data"] += size
        elif sec.startswith(".bss"): tot["bss"] += size
    def word(s):
        lines = (out / f"v{v}.{s}.txt").read_text().splitlines()
        hexes = [l.split()[1] for l in lines if re.match(r"^\s*0000 ", l)]
        b = bytes.fromhex(hexes[0])
        return int.from_bytes(b[:4], "little")
    su = [l.rsplit("\t", 2) for l in (out / f"v{v}.su").read_text().splitlines()]
    fn, frame, _ = max(su, key=lambda x: int(x[1]))
    fn = re.sub(r"^.*?:\d+:\d+:", "", fn)
    fn = re.sub(r" \[with .*", "", fn)[:90]
    wr.writerow([n, tot["code"], tot["rodata"], tot["data"], tot["bss"], word("fusion"), word("mekf"), frame, fn])
PY
