"""Pinned SDK FCR writer inventory and lazy-context correspondence.

This is a conditional preservation audit. The architecture leaves FCR's
reset value undefined, and the task's first use inherits the active value.
Neither a zeroed task stack nor absence of fesetround proves initialization.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import subprocess

from tools.stability.ou3_alt_contraction import target_scalar_profile as S

QUALIFICATION='OU3_ALT_PINNED_FCR_CONTEXT_AUDIT_V1'
EXPECTED_WRITERS={
    ('lib/libesp_gdbstub.a','esp_gdbstub_set_register'),
    ('lib/libxtensa.a','_xt_coproc_exc'),
    *(('lib/libxt_hal.a',name) for name in
      ('xthal_restore_cp0','xthal_restore_cp0_nw',
       'xthal_restore_cpregs','xthal_restore_cpregs_nw')),
}


def instructions(disassembly):
    """Keep archive member, symbol, address and decoded instruction."""
    member=symbol=''
    result=[]
    for line in disassembly.splitlines():
        if ':     file format' in line: member=line.split(':',1)[0]
        match=re.match(r'[0-9a-f]+ <([^>]+)>:',line)
        if match: symbol=match[1]
        match=re.match(r'\s*([0-9a-f]+):\s+([0-9a-f]+)\s+(.+)',line)
        if match:
            result.append((member,symbol,int(match[1],16),match[2],
                           re.sub(r'\s+',' ',match[3].strip())))
    return result


def context_sites(code):
    ins=instructions(code)
    sites={address:(encoding,op) for _,symbol,address,encoding,op in ins
           if symbol=='_xt_coproc_exc'}
    expected={0x207:'e33e80',0x20a:'0239', # FCR -> old save-area word zero
              0x2c7:'011f32',0x2cd:'028307',0x2d0:'003906', # full-saved gate
              0x2e7:'0238',0x2e9:'f3e830', # new save-area word zero -> FCR
              0x3b8:'021f22',0x3bb:'e50207'} # first use skips restore
    if any(sites.get(address,(None,))[0]!=encoding
           for address,encoding in expected.items()):
        raise ValueError('pinned lazy FCR context operation changed')
    if [x[4] for x in ins if x[1]=='_xt_coproc_exc' and 'wur.fcr' in x[4]]!=['wur.fcr a3']:
        raise ValueError('additional FCR writer in lazy coprocessor handler')
    return {hex(address):sites[address][1] for address in expected}


def audit(toolchain:Path,sdk:Path,rom:Path):
    objdump=toolchain/'bin/xtensa-esp32s3-elf-objdump'
    def scan(path):
        code=S.run([objdump,'-d',path])
        ins=instructions(code)
        return {'path':str(path.relative_to(sdk)), 'sha256':S.sha(path),
                'FCR_accesses':[{'member':m,'symbol':s,'address':hex(a),
                                'encoding':e,'operation':op}
                               for m,s,a,e,op in ins if '.fcr' in op]}
    archives=sorted((sdk/'lib').glob('*.a'))+sorted((sdk/'qio_qspi').glob('*.a'))
    with ThreadPoolExecutor(max_workers=8) as pool: inventory=list(pool.map(scan,archives))
    writers={(item['path'],access['symbol']) for item in inventory
             for access in item['FCR_accesses'] if access['operation'].startswith('wur.fcr')}
    if writers!=EXPECTED_WRITERS:
        raise ValueError('pinned SDK FCR writer inventory changed')
    sites=context_sites(S.run([objdump,'-d',sdk/'lib/libxtensa.a']))
    startup=instructions(S.run([objdump,'-d',sdk/'qio_qspi/libfreertos.a']))
    task={address:encoding for _,symbol,address,encoding,_ in startup
          if symbol=='pxPortInitialiseStack'}
    if any(task.get(a)!=e for a,e in {0x3:'345020',0x6:'059c',
                                    0x36:'0759',0x38:'1759',0x3a:'2799'}.items()):
        raise ValueError('task coprocessor flags initialization changed')
    rom_access=[x for x in instructions(S.run([objdump,'-d',rom])) if '.fcr' in x[4]]
    if rom_access: raise ValueError('ROM FCR access requires a new startup analysis')
    return {'qualification':QUALIFICATION,'isa_url':S.ISA,
        'isa_reset_table':'Table 194, FCR user register 232 (printed page 317)',
        'architectural_FCR_reset_value':'undefined',
        'sdk_memory_variant':'qio_qspi','sdk_archive_count':len(inventory),
        'compiler_objdump_sha256':S.sha(objdump.resolve()),
        'sdk_archives':inventory,'rom_elf_sha256':S.sha(rom),
        'rom_FCR_access_count':len(rom_access),'lazy_context_sites':sites,
        'task_first_FPU_use_inherits_active_FCR':True,
        'saved_context_restores_its_original_FCR':True,
        'RNE_preserved_if_all_active_and_saved_contexts_are_RNE':True,
        'required_execution_premises':[
            'both cores enter the admitted execution with FCR.RM=0',
            'saved coprocessor areas are valid and unmodified',
            'no debugger or other code changes FCR.RM'],
        'rounding_mode_runtime_initialization_and_preservation_proved':False,
        'whole_firmware_link_and_callgraph_qualified':False,
        'source_reference_urls':[
          'https://github.com/espressif/esp-idf/blob/v5.5.2/components/freertos/FreeRTOS-Kernel/portable/xtensa/port.c',
          'https://github.com/espressif/esp-idf/blob/v5.5.2/components/xtensa/xtensa_vectors.S'],
        'storage_search_allowed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('toolchain','sdk','rom','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args(); result=audit(a.toolchain,a.sdk,a.rom)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('sdk_archive_count',
        'task_first_FPU_use_inherits_active_FCR',
        'rounding_mode_runtime_initialization_and_preservation_proved')}))


if __name__=='__main__': main()
