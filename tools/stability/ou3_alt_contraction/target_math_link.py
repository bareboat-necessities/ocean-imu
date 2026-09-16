"""Inspect pinned exp/log/sqrt namespace resolution with the SDK ROM script.

This is a static standalone link probe, not an executable Arduino firmware.
It resolves real library objects without supplying replacement math functions.
The SDK flash instruction segment contains every text/literal section; the
SDK ROM divide entry remains unchanged. Bare libc's __getreent diagnostic is
retained in the report, so success cannot imply complete runtime linkage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

from tools.stability.ou3_alt_contraction import target_qaxis_exp as Q
from tools.stability.ou3_alt_contraction import target_scalar_profile as S

QUALIFICATION = 'OU3_ALT_PINNED_MATH_NAMESPACE_LINK_V1'
PROBE = '''#include <cmath>
extern "C" float alt_math_probe(float e,float l,float s) {
    return std::exp(e)+std::log(l)+std::sqrt(s);
}
'''
LAYOUT = '''SECTIONS
{
  . = 0x42000020;
  .text : { *(.literal .literal.*) *(.text .text.*) }
  .rodata : { *(.rodata .rodata.*) }
  .data : { *(.data .data.*) }
  .bss : { *(.bss .bss.*) *(COMMON) }
}
'''


def _run(args):
    completed = subprocess.run([str(a) for a in args], text=True,
                               capture_output=True)
    if completed.returncode:
        raise RuntimeError('target math namespace probe failed: '+completed.stderr)
    return completed


def audit(toolchain: Path, sdk: Path):
    compiler = toolchain/'bin/xtensa-esp32s3-elf-g++'
    nm = toolchain/'bin/xtensa-esp32s3-elf-nm'
    objdump = toolchain/'bin/xtensa-esp32s3-elf-objdump'
    libm = Path(_run([compiler, '-print-file-name=libm.a']).stdout.strip())
    if S.sha(libm) != Q.LIBM_SHA256:
        raise ValueError('unqualified target libm bytes')
    gcc = Path(_run([compiler, '-print-libgcc-file-name']).stdout.strip())
    rom_ld = sdk/'ld/esp32s3.rom.libgcc.ld'
    if 'org = 0x42000020' not in (sdk/'ld/memory.ld').read_text():
        raise ValueError('SDK flash instruction segment changed')
    # No SDK script redirects the actual exp/log/sqrt symbols being checked.
    symbols = ('expf', 'logf', 'sqrtf', '__ieee754_expf', '__ieee754_logf', '__ieee754_sqrtf')
    for ld in (sdk/'ld').glob('*.ld'):
        if any(re.search(r'\b'+re.escape(name)+r'\s*=', ld.read_text()) for name in symbols):
            raise ValueError('SDK math-symbol override requires qualification')
    flags = shlex.split((sdk/'flags/cpp_flags').read_text())
    flags += ['-Os', '-funroll-loops', '-fno-finite-math-only']
    with tempfile.TemporaryDirectory(prefix='ou3-math-link-') as directory:
        tmp = Path(directory)
        source, layout, obj, elf, link_map = [tmp/name for name in
            ('probe.cpp', 'probe.ld', 'probe.o', 'probe.elf', 'probe.map')]
        source.write_text(PROBE)
        layout.write_text(LAYOUT)
        _run([compiler, *flags, '-ffile-prefix-map='+str(tmp)+'=/ou3-math-link',
              '-c', source, '-o', obj])
        link = _run([compiler, '-nostartfiles', '-Wl,-e,alt_math_probe',
                     '-Wl,--gc-sections', '-Wl,-T,'+str(layout),
                     '-Wl,-T,'+str(rom_ld), '-Wl,-Map,'+str(link_map), obj,
                     '-o', elf, '-lm', '-lc', '-lgcc'])
        names = {}
        for line in _run([nm, '-n', elf]).stdout.splitlines():
            match = re.match(r'^([0-9a-f]+)\s+(\w)\s+(\S+)$', line)
            if match and match[3] in symbols+('__divsf3',):
                names[match[3]] = {'address': '0x'+match[1], 'kind': match[2]}
        if names.get('__divsf3') != {'address':'0x40002274', 'kind':'A'}:
            raise ValueError('namespace probe does not use the SDK ROM divide entry')
        for name in symbols:
            if name not in names or names[name]['kind'] != 'T':
                raise ValueError('math namespace did not resolve '+name)
            if not 0x42000020 <= int(names[name]['address'], 16) < 0x42800000:
                raise ValueError('math code placed outside SDK instruction segment')
        mapping = link_map.read_text()
        for member in ('libm_a-wf_exp.o', 'libm_a-ef_exp.o', 'libm_a-wf_log.o',
                       'libm_a-ef_log.o', 'libm_a-wf_sqrt.o', '_sqrtf.o'):
            if '('+member+')' not in mapping:
                raise ValueError('namespace supplier missing: '+member)
        dis = _run([objdump, '-d', elf]).stdout
        if [r[1] for r in S.body(dis, '__ieee754_sqrtf')] != S.SQRT_ENCODINGS:
            raise ValueError('linked sqrt differs from audited IEEE sequence')
        # Preserve diagnostics without machine-specific installation paths.
        diagnostic = link.stderr.replace(str(toolchain), '/toolchain')
        elf_sha = S.sha(elf)
    return {
        'qualification': QUALIFICATION,
        'compiler_driver_sha256': S.sha(compiler.resolve()),
        'sdk_cpp_flags_sha256': S.sha(sdk/'flags/cpp_flags'),
        'sdk_rom_linker_sha256': S.sha(rom_ld),
        'libm_sha256': S.sha(libm), 'libgcc_sha256': S.sha(gcc),
        'probe_source_sha256': hashlib.sha256(PROBE.encode()).hexdigest(),
        'probe_layout_sha256': hashlib.sha256(LAYOUT.encode()).hexdigest(),
        'standalone_probe_elf_sha256': elf_sha,
        'resolved_symbols': names,
        'standalone_exp_log_sqrt_namespace_resolution_qualified': True,
        'linked_sqrt_matches_audited_IEEE_sequence': True,
        'division_resolves_to_audited_SDK_ROM_entry': True,
        'SDK_link_scripts_override_exp_log_sqrt': False,
        'non_math_runtime_linker_diagnostic': diagnostic.strip(),
        'standalone_probe_runtime_complete': False,
        'whole_firmware_link_resolution_qualified': False,
        'target_execution_observed': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--toolchain', type=Path, required=True)
    parser.add_argument('--sdk', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.toolchain.resolve(), args.sdk.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result['resolved_symbols'], sort_keys=True))


if __name__ == '__main__':
    main()
