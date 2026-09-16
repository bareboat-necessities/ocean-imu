"""Audit the pinned Xtensa scalar operations and the actual ROM divide entry.

The ISA supplies operation semantics; matching the ROM/program bytes supplies
the implementation identity. FCR.RM=0 is an explicit execution premise, not an
inference from a host test. This does not qualify arbitrary compiler reductions.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

QUALIFICATION='OU3_ALT_XTENSA_SCALAR_PROFILE_V1'
ISA='https://www.cadence.com/content/dam/cadence-www/global/en_US/documents/tools/silicon-solutions/compute-ip/isa-summary.pdf'
ROM_URL='https://github.com/espressif/esp-rom-elfs/releases/download/20241011/esp-rom-elfs-20241011.tar.gz'
ROM_ARCHIVE_SHA='921f000164a421c7628fbfee55b173384aafaa51883adc65cd27bf9b0af9e9a9'
# These instruction encodings are taken from the audited target program, not
# from the ISA prose. Both sequences use the final IEEE correction operation.
DIV_ENCODINGS='002136 fa1250 fa2350 fa3270 fa42b0 fa5130 6a5430 fa6300 fa7200 fa21b0 6a6560 fa5130 fa0030 fa8260 6a5460 6a0830 fa71d0 6a6560 6a8400 fa3130 6a3460 6a0860 fa2260 6a6360 6a2400 fa07f0 fa67e0 7a0260 fa2040 f01d'.split()
SQRT_ENCODINGS='002136 fa1250 fa2190 fa3030 6a3220 fa41b0 fa0330 fa40e0 6a0340 fa31b0 fa5360 6a2020 fa0030 fa6030 fa7030 6a0520 6a6240 fa4330 6a7420 6a3000 6a4620 fa2760 6a0320 6a7470 fa21c0 fa11b0 6a1000 fa3760 fa02f0 fa32e0 7a0130 fa2040 f01d'.split()


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args):
    return subprocess.run([str(x) for x in args],check=True,text=True,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout


def body(disassembly,symbol):
    match=re.search(r'^([0-9a-f]+) <'+re.escape(symbol)+r'>:\n',disassembly,re.M)
    if not match: raise ValueError('missing target symbol '+symbol)
    result=[]
    for line in disassembly[match.end():].splitlines():
        item=re.match(r'\s*([0-9a-f]+):\s+([0-9a-f]+)\s+(.+)',line)
        if not item: break
        result.append((int(item[1],16),item[2],item[3].strip()))
        if item[3].strip()=='retw.n': break
    return result


def audit(toolchain,sdk,rom):
    compiler=toolchain/'bin/xtensa-esp32s3-elf-g++'
    objdump=toolchain/'bin/xtensa-esp32s3-elf-objdump'
    gcc=Path(run([compiler,'-print-libgcc-file-name']).strip())
    archive=run([objdump,'-dr',gcc])
    rom_code=run([objdump,'-d',rom])
    rd=body(rom_code,'__divsf3'); gd=body(archive,'__divsf3')
    sr=body(archive,'__ieee754_sqrtf')
    if [x[1] for x in rd]!=DIV_ENCODINGS or [x[1] for x in gd]!=DIV_ENCODINGS:
        raise ValueError('ROM/archive divide differs from audited IEEE sequence')
    if [x[1] for x in sr]!=SQRT_ENCODINGS:
        raise ValueError('archive square root differs from audited IEEE sequence')
    ld=sdk/'ld/esp32s3.rom.libgcc.ld'
    if not re.search(r'^__divsf3\s*=\s*0x40002274;',ld.read_text(),re.M):
        raise ValueError('SDK ROM divide entry changed')
    thunk=run([objdump,'-d','--start-address=0x40002274','--stop-address=0x4000227a',rom])
    if not (rd[0][0]==0x40056124 and
            re.search(r'l32r\s+a9,.*\(40056124 <__divsf3>\)',thunk) and
            re.search(r'jx\s+a9',thunk)):
        raise ValueError('SDK-selected ROM thunk no longer reaches audited divide')
    # The SDK does not replace the sqrt implementation by an unaudited ROM
    # address. The final firmware link remains a separate checked obligation.
    if any('__ieee754_sqrtf' in p.read_text() for p in (sdk/'ld').glob('*.ld')):
        raise ValueError('new SDK square-root symbol override requires review')
    return {
        'qualification':QUALIFICATION,
        'isa_url':ISA,'isa_sections':['4.3.11.2','4.3.11.3','4.3.11.5','8.3.159','8.3.161'],
        'rom_archive_url':ROM_URL,'rom_archive_sha256':ROM_ARCHIVE_SHA,
        'rom_elf_sha256':sha(rom),'libgcc_sha256':sha(gcc),'sdk_rom_linker_sha256':sha(ld),
        'rom_divide_entry':'0x40002274','rom_divide_body':'0x40056124',
        'ROM_thunk_and_IEEE_divide_sequence_correspondence_closed':True,
        'archive_IEEE_sqrt_sequence_correspondence_closed':True,
        'scalar_add_sub_mul_and_fused_madd_msub_IEEE_semantics_closed':True,
        'gradual_underflow_IEEE_semantics_closed':True,
        'rounding_premise':'FCR.RM=0 throughout the admitted execution',
        'rounding_mode_runtime_initialization_and_preservation_proved':False,
        'division_correctly_rounded_under_rounding_premise':True,
        'sqrt_correctly_rounded_under_rounding_and_link_premises':True,
        'whole_firmware_link_resolution_qualified':False,
        'all_expression_contractions_and_reductions_qualified':False,
        'target_execution_observed':False,'storage_search_allowed':False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--toolchain',type=Path,required=True)
    p.add_argument('--sdk',type=Path,required=True)
    p.add_argument('--rom',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); report=audit(a.toolchain,a.sdk,a.rom)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('qualification','ROM_thunk_and_IEEE_divide_sequence_correspondence_closed','storage_search_allowed')}))


if __name__=='__main__': main()
