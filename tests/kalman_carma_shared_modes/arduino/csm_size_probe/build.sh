#!/bin/bash
# ESP32-S3 size / stack probe of the standalone shared-modes QMEKF (see
# csm_size_probe.cpp).  Usage: build.sh <xtensa toolchain bin dir> [out dir]
# The toolchain is Espressif's crosstool-NG release (esp-14.2.0_20241119, the
# GCC of arduino-esp32 3.3.x); flags follow arduino-esp32's C++ defaults.
# Prints CSV: sizes, the largest single stack frame, and the worst static
# call-chain stack from csm_probe_step (-fcallgraph-info=su; calls through
# function pointers and the libraries' own frames are not counted, so it is
# a lower bound of the real high-water mark, not a measurement).
set -euo pipefail
TC=${1:?xtensa toolchain bin dir}
OUT=${2:-$(mktemp -d)}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../.." && pwd)
EIGEN=${EIGEN_DIR:-$ROOT/third_party/eigen}
[ -f "$EIGEN/Eigen/Dense" ] || EIGEN=/usr/include/eigen3
CXX="$TC/xtensa-esp32s3-elf-g++"; OBJDUMP="$TC/xtensa-esp32s3-elf-objdump"
for v in 0 1; do
  (cd "$OUT" && "$CXX" -std=gnu++2a -Os -mlongcalls -ffunction-sections -fdata-sections \
      -fexceptions -fno-rtti -fstack-usage -fcallgraph-info=su -DCSM_PROBE_VARIANT=$v \
      -I"$ROOT/src" -isystem "$EIGEN" -c "$HERE/csm_size_probe.cpp" -o "$OUT/v$v.o")
  "$OBJDUMP" -h "$OUT/v$v.o" > "$OUT/v$v.h.txt"
  for s in top wrapper mekf; do
    "$OBJDUMP" -s -j ".rodata.csm_probe_sizeof_$s" "$OUT/v$v.o" > "$OUT/v$v.$s.txt"
  done
done
python3 - "$OUT" <<'PY'
import csv, re, sys
from collections import defaultdict
from pathlib import Path
out = Path(sys.argv[1])
wr = csv.writer(sys.stdout, lineterminator="\n")
wr.writerow(["variant", "code_bytes", "rodata_bytes", "sizeof_top", "sizeof_wrapper", "sizeof_qmekf",
             "max_stack_frame", "max_stack_function", "update_chain_stack", "update_chain_depth"])
for v, n in enumerate(["A41", "CSM"]):
    tot = {"code": 0, "rodata": 0}
    for line in (out / f"v{v}.h.txt").read_text().splitlines():
        f = line.split()
        if len(f) < 3 or not f[0].isdigit():
            continue
        sec, size = f[1], int(f[2], 16)
        if sec.startswith((".text", ".literal")): tot["code"] += size
        elif sec.startswith(".rodata"): tot["rodata"] += size
    def word(s):
        lines = (out / f"v{v}.{s}.txt").read_text().splitlines()
        hexes = [l.split()[1] for l in lines if re.match(r"^\s*0000 ", l)]
        return int.from_bytes(bytes.fromhex(hexes[0])[:4], "little")
    su = [l.rsplit("\t", 2) for l in (out / f"v{v}.su").read_text().splitlines()]
    fn, frame, _ = max(su, key=lambda x: int(x[1]))
    fn = re.sub(r" \[with .*", "", re.sub(r"^.*?:\d+:\d+:", "", fn))[:80]
    # call graph (VCG): node: { title: "x" label: "name\n...\nN bytes (static)" } edge: { sourcename: "a" targetname: "b" }
    ci = (out / f"v{v}.ci").read_text(errors="replace")
    stack, edges = {}, defaultdict(set)
    for m in re.finditer(r'node: \{ title: "([^"]+)" label: "([^"]*)"', ci):
        b = re.search(r"(\d+) bytes", m.group(2))
        stack[m.group(1)] = int(b.group(1)) if b else 0
    for m in re.finditer(r'edge: \{ sourcename: "([^"]+)" targetname: "([^"]+)"', ci):
        edges[m.group(1)].add(m.group(2))
    root = next((t for t in stack if "csm_probe_step" in t), None)
    memo, onpath = {}, set()
    def worst(t):
        if t in memo: return memo[t]
        if t in onpath: return (0, 0)
        onpath.add(t)
        best = (0, 0)
        for c in edges.get(t, ()):
            w = worst(c)
            if w[0] > best[0]: best = w
        onpath.discard(t)
        memo[t] = (stack.get(t, 0) + best[0], best[1] + 1)
        return memo[t]
    chain, depth = worst(root) if root else (float("nan"), 0)
    wr.writerow([n, tot["code"], tot["rodata"], word("top"), word("wrapper"), word("mekf"), frame, fn, chain, depth])
PY
