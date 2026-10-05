"""Parity-resolved homogeneous error-factor product for the planar 20-s word.

FINITE DIAGNOSTIC ONLY. Reports full, even and odd coordinate gains so an
expanding unobservable/integrator direction is not hidden inside a Mahony-only
surrogate contraction.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .planar_compatibility_quotient_mean import word,SCOPE
from .planar_parity import EVEN,ODD

def diagnostic(path,start=40000,end=44000):
 M,_,_=word(path,start,end)
 def sv(idx):return np.linalg.svd(M[np.ix_(idx,idx)],compute_uv=False)
 sf=np.linalg.svd(M,compute_uv=False);se=sv(EVEN);so=sv(ODD)
 return {"qualification":"OU3_PLANAR_PARITY_MEAN_TANGENT_V1","result_type":"FINITE DIAGNOSTIC ONLY",
  "word_samples":[start,end],"full_gain_2norm":float(sf[0]),
  "operator_scope":SCOPE,"complete_nonlinear_mean_derivative_computed":False,
  "even_gain_2norm":float(se[0]),"odd_gain_2norm":float(so[0]),
  "even_leading_singular_values":se[:8].tolist(),"odd_leading_singular_values":so[:8].tolist(),
  "mahony_only_gain_is_complete_mean_gain":False,
  "complete_mean_contraction_verified":False,
  "finite_homogeneous_euclidean_contraction":bool(max(se[0],so[0])<1),
  "euclidean_gain_proves_instability":False,
  "all_time_magnetic_service_verified":False,"theorem_closed":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=diagnostic(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
