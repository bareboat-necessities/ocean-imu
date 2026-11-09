#!/usr/bin/env python3
"""Compile one deterministic trace harness against pristine main and candidate.

Usage: --baseline <checkout at 1fe4711> --eigen <Eigen include root>.
The source, random stream, physical inputs, startup, adaptation, temperature and
magnetic outage are identical. Timings are excluded. Existing gate limits are
unchanged; this additional float trace comparison uses 1e-6 + 2e-6*abs(reference).
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3'));ap.add_argument('--output',type=Path,default=ROOT/'reports/results/tfg_shared_modes/refactor_regression.json');args=ap.parse_args()
source=ROOT/'tests/kalman_carma_shared_modes/refactor_trace.cpp'
result=dict(baseline_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.baseline,text=True).strip(),
 harness_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),absolute_tolerance=1e-6,relative_tolerance=2e-6,filters={})
with tempfile.TemporaryDirectory(prefix='tfg-refactor-') as temp:
 work=Path(temp)
 for kind,name in enumerate(['QMEKF-CARMA','TFG-CARMA','QMEKF-shared','OU-III']):
  arrays=[]
  for variant,root in [('baseline',args.baseline),('candidate',ROOT)]:
   binary=work/f'{kind}_{variant}'
   cmd=['g++','-std=c++20','-O1','-DEIGEN_UNROLLING_LIMIT=0',f'-DREFACTOR_FILTER={kind}','-I',str(root/'src'),'-isystem',str(args.eigen),str(source),'-o',str(binary)]
   subprocess.run(cmd,check=True)
   trace=work/f'{kind}_{variant}.csv'
   with trace.open('w') as out:subprocess.run([str(binary)],stdout=out,check=True)
   arrays.append(np.loadtxt(trace,delimiter=','))
  a,b=arrays
  if not np.array_equal(a[:,:2],b[:,:2]):raise AssertionError(f'{name}: startup chronology changed')
  nanmask=np.isnan(a)==np.isnan(b)
  delta=np.abs(a-b);scaled=delta/(1e-6+2e-6*np.abs(a));peak=float(np.nanmax(scaled));ok=bool(nanmask.all() and peak<=1.)
  result['filters'][name]=dict(samples=len(a),identical=bool(np.array_equal(a,b,equal_nan=True)),max_absolute_difference=float(np.nanmax(delta)),max_scaled_difference=peak,passed=ok)
  print(name,result['filters'][name],flush=True)
  if not ok:raise AssertionError(name)
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
