#!/usr/bin/env python3
"""Build matched AtomS3R applications and the on-device measurement harness.

Install esp32:esp32@3.3.7, M5Unified@0.2.13, M5GFX@0.2.19 and Eigen@0.3.2
first. This records linked sizes; serial timing/heap/stack requires hardware.
"""
import argparse,json,re,subprocess,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--arduino-cli',default='arduino-cli')
 ap.add_argument('--config-file',type=Path)
 ap.add_argument('--build-root',type=Path,required=True)
 ap.add_argument('--results',type=Path,default=ROOT/'reports/results/tfg_shared_modes')
 args=ap.parse_args();args.results.mkdir(parents=True,exist_ok=True)
 cli=[args.arduino_cli]
 if args.config_file:cli+=['--config-file',str(args.config_file)]
 def inventory(command):
  r=subprocess.run(cli+command,capture_output=True,text=True,check=True)
  return json.loads(r.stdout)
 cores={x['id']:x.get('installed_version') for x in inventory(['core','list','--format','json'])['platforms']}
 libs={x['library']['name']:x['library']['version'] for x in inventory(['lib','list','--format','json'])['installed_libraries']}
 if cores.get('esp32:esp32')!='3.3.7':raise SystemExit('Install esp32:esp32@3.3.7 for the matched measurement.')
 for name,version in {'M5Unified':'0.2.13','M5GFX':'0.2.19','Eigen':'0.3.2','ArxTypeTraits':'0.3.2'}.items():
  if libs.get(name)!=version:raise SystemExit('Install '+name+'@'+version+' for the matched measurement.')
 record=dict(platform='esp32:esp32@3.3.7',on_device_timing=False,
             cli=inventory(['version','--format','json']),installed_cores=cores,
             installed_libraries=libs,builds=[])
 variants=[('qmekf','sensors/full_marine_ins/atomS3R_ins_carma_shared_modes',''),
           ('tfg','sensors/full_marine_ins/atomS3R_ins_carma_shared_modes','-DOCEAN_IMU_TFG_SHARED_MODES=1'),
           ('probe','tests/kalman_carma_shared_modes/arduino/tfg_shared_runtime_probe','-DPROBE_VARIANT=3')]
 for name,sketch,extra in variants:
  dest=args.build_root.resolve()/name;dest.mkdir(parents=True,exist_ok=True)
  command=cli+['compile','--fqbn','esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc',
   '--jobs','1','--library',str(ROOT),'--build-path',str(dest),'--output-dir',str(dest/'output'),
   '--build-property','compiler.cpp.extra_flags=-g0 -funroll-loops -fno-finite-math-only '+extra,str(ROOT/sketch)]
  start=time.monotonic();log=dest/'build.log'
  with log.open('w') as stream:r=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
  text=log.read_text();item=dict(variant=name,exit_code=r.returncode,command=command,
                              elapsed_seconds=time.monotonic()-start,log_tail=text[-1800:])
  for pattern,key in [(r'Sketch uses (\d+) bytes','flash_bytes'),(r'Global variables use (\d+) bytes','static_ram_bytes')]:
   match=re.search(pattern,text)
   if match:item[key]=int(match[1])
  record['builds'].append(item)
  (args.results/'firmware_builds.json').write_text(json.dumps(record,indent=2)+'\n')
  print(name,'exit',r.returncode,'flash',item.get('flash_bytes'),'RAM',item.get('static_ram_bytes'),flush=True)
  if r.returncode:raise SystemExit(r.returncode)

if __name__=='__main__':main()
