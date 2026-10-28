#!/bin/bash
# ESP32-S3 link and call-graph check of the size probe for OU-III (0), CARMA3
# C41 (1) and TFG-CARMA3 (7).  Usage: callgraph.sh <xtensa toolchain bin dir> [out dir]
# Same toolchain and flags as build.sh, plus -fcallgraph-info=su.  Writes
# esp32s3_link.csv and esp32s3_callgraph.csv to the out dir:
#   * the link of each probe object with link_main.cpp against the toolchain's
#     newlib/libstdc++ (default linker script, --gc-sections): undefined
#     references, warnings other than newlib's OS stubs (which ESP-IDF
#     supplies on target), and the linked text/data/bss;
#   * per sample-path entry and filter method, the deepest static call chain
#     and the reachability of heap allocation, soft-float double arithmetic
#     and C++ throw (callgraph.py).
# Build-time evidence only: no on-target timing.
set -euo pipefail
TC=${1:?xtensa toolchain bin dir}
OUT=${2:-$(mktemp -d)}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../.." && pwd)
EIGEN=${EIGEN_DIR:-/usr/include/eigen3}
CXX="$TC/xtensa-esp32s3-elf-g++"
FLAGS=(-std=gnu++2a -Os -mlongcalls -ffunction-sections -fdata-sections -fexceptions -fno-rtti)
(cd "$OUT" && "$CXX" "${FLAGS[@]}" -c "$HERE/link_main.cpp" -o "$OUT/link_main.o")
for v in 0 1 7; do
  (cd "$OUT" && "$CXX" "${FLAGS[@]}" -fstack-usage -fcallgraph-info=su -DCARMA3_PROBE_VARIANT=$v \
      -I"$ROOT/src" -isystem "$EIGEN" -c "$HERE/carma3_size_probe.cpp" -o "$OUT/v$v.o")
  rc=0
  "$CXX" -mlongcalls -Wl,--gc-sections -Wl,--warn-unresolved-symbols "$OUT/link_main.o" "$OUT/v$v.o" \
      -o "$OUT/v$v.elf" > "$OUT/v$v.link.log" 2>&1 || rc=$?
  echo "$rc" > "$OUT/v$v.link.rc"
  "$TC/xtensa-esp32s3-elf-size" -B "$OUT/v$v.elf" | tail -1 > "$OUT/v$v.size.txt" || true
done
names=(OU3 A41 x x x x x T41)
classes=(Kalman3D_Wave_OU_III Kalman3D_Wave_CARMA3 x x x x x Kalman3D_Wave_TFG_CARMA3)
{
echo "variant,link_exit,undefined_references,other_link_warnings,text_bytes,data_bytes,bss_bytes"
for v in 0 1 7; do
  log="$OUT/v$v.link.log"
  undef=$(grep -c "undefined reference" "$log" || true)
  other=$(grep -cvE "is not implemented and will always fail|does not take linker garbage collection" "$log" || true)
  read -r text data bss _ < "$OUT/v$v.size.txt"
  echo "${names[$v]},$(cat "$OUT/v$v.link.rc"),$undef,$other,$text,$data,$bss"
done
} > "$OUT/esp32s3_link.csv"
{
echo "variant,root,worst_stack_chain_bytes,chain_frames,reaches_heap_allocation,reaches_double_arithmetic,reaches_throw,largest_frame_in_chain"
for v in 0 1 7; do
  python3 "$HERE/callgraph.py" "${names[$v]}" "$OUT/v$v.ci" "${classes[$v]}"
done
} > "$OUT/esp32s3_callgraph.csv"
echo "wrote $OUT/esp32s3_link.csv and $OUT/esp32s3_callgraph.csv"
