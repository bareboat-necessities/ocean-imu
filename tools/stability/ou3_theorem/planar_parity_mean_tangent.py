"""Parity-resolved point mean tangent for the planar 20-s word.

FINITE DIAGNOSTIC ONLY. Reports full, even and odd coordinate gains so an
expanding unobservable/integrator direction is not hidden inside a Mahony-only
surrogate contraction.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import reset_matrix
from .planar_parity import EVEN,ODD

def diagnostic(path,start=40000,end=44000):
 M=np.eye(21)
 for kind,k,a in records(path):
  if k<=start:continue
  if k>end:break
  if kind==1:M=expand(a[225:450])@M
  elif kind==9:
   H=a[225:288].reshape(3,21);K=a[297:360].reshape(21,3);M=(np.eye(21)-K@H)@M
  elif kind==5:M=reset_matrix(a[225:228])@M
 def sv(idx):return np.linalg.svd(M[np.ix_(idx,idx)],compute_uv=False)
 sf=np.linalg.svd(M,compute_uv=False);se=sv(EVEN);so=sv(ODD)
 return {"qualification":"OU3_PLANAR_PARITY_MEAN_TANGENT_V1","result_type":"FINITE DIAGNOSTIC ONLY",
  "word_samples":[start,end],"full_gain_2norm":float(sf[0]),
  "even_gain_2norm":float(se[0]),"odd_gain_2norm":float(so[0]),
  "even_leading_singular_values":se[:8].tolist(),"odd_leading_singular_values":so[:8].tolist(),
  "mahony_only_gain_is_complete_mean_gain":False,
  "complete_mean_contraction_verified":bool(max(se[0],so[0])<1),
  "quotient_or_transverse_projection_required":bool(max(se[0],so[0])>=1),
  "all_time_magnetic_service_verified":False,"theorem_closed":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=diagnostic(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
