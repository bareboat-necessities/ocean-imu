#!/usr/bin/env bash
# Usage: bash tools/marine_mag_resources.sh TOOLCHAIN_BIN SCRATCH_DIR
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
bin=$(realpath "$1")
out=$(realpath -m "$2")
mkdir -p "$out"
cd "$out"
src="$root/tests/kalman_tfg_shared_modes/magnetic_study"
eigen=${EIGEN_DIR:-$root/third_party/eigen}
[[ -f "$eigen/Eigen/Dense" ]] || eigen=/usr/include/eigen3
flags=(-std=c++20 -I"$root/src" -isystem "$eigen")
for unit in resources SelectiveMag; do
    "$bin/xtensa-esp32s3-elf-g++" "${flags[@]}" -Os -mlongcalls -ffunction-sections -fdata-sections \
        -fstack-usage -c "$src/$unit.cpp" -o "$unit.o"
done
"$bin/xtensa-esp32s3-elf-g++" --version | head -1
"$bin/xtensa-esp32s3-elf-objdump" -s -j .rodata.learner_bytes -j .rodata.gate_bytes resources.o
"$bin/xtensa-esp32s3-elf-size" resources.o SelectiveMag.o
awk -F '\t' 'NF>=3 && $2+0 > max {max=$2+0; line=$0} END {print line}' resources.su
awk -F '\t' 'NF>=3 && $2+0 > max {max=$2+0; line=$0} END {print line}' SelectiveMag.su
g++ "${flags[@]}" -O3 "$src/resources.cpp" "$src/SelectiveMag.cpp" -o timing
./timing
"$bin/xtensa-esp32s3-elf-g++" "${flags[@]}" -Os -mlongcalls -ffunction-sections -fdata-sections \
    -fstack-usage -c "$src/packet-resources.cpp" -o packet-resources.o
"$bin/xtensa-esp32s3-elf-size" packet-resources.o
awk -F '\t' 'NF>=3 && $2+0 > max {max=$2+0; line=$0} END {print line}' packet-resources.su
g++ "${flags[@]}" -O3 "$src/packet-resources.cpp" -o packet-timing
./packet-timing
