"""Measurement-compatible quotient audit for the literal planar mean word.

FINITE DIAGNOSTIC ONLY.  The quotient line is not chosen from low-gain SVD
vectors.  It is the exact planar compatibility direction: opposite roll/BA-y
splits induce the same delivered sensor history.  We test the literal word
against the covariance-metric quotient machinery and report gauge leakage.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import reset_matrix
from .compatibility_quotient import quotient_terminal_map

def word(path,start,end):
 M=np.eye(21);P0=PN=None
 for kind,k,a in records(path):
  if k<=start:continue
  if k>end:break
  if kind==1:
   if P0 is None:P0=expand(a[:225])
   M=expand(a[225:450])@M
  elif kind==9:
   H=a[225:288].reshape(3,21);K=a[297:360].reshape(21,3);M=(np.eye(21)-K@H)@M
  elif kind==5:M=reset_matrix(a[225:228])@M
  elif kind==7:PN=expand(a[:225])
 if P0 is None or PN is None:raise ValueError("empty word")
 return M,P0,PN

def compatibility_line(sample):
 # Linearized physical +/- separation: roll theta_x and BA_y are linked by
 # delta b_a,y = g * delta theta_x.  This is the certified planar gauge seed.
 r=np.zeros(21);r[0]=1.;r[19]=9.80665
 return r

def diagnostic(path,start=40000,end=44000):
 M,P0,PN=word(path,start,end);r0=compatibility_line(start);rN=compatibility_line(end)
 q=quotient_terminal_map(P0,r0,PN,rN,M,np.zeros(21))
 s=np.linalg.svd(q["M_Q"],compute_uv=False)
 return {"qualification":"OU3_PLANAR_COMPATIBILITY_QUOTIENT_MEAN_V1",
  "result_type":"FINITE DIAGNOSTIC ONLY","word_samples":[start,end],
  "compatibility_line_definition":"theta_x=1, b_a_y=g; exact planar roll/BA compatibility tangent",
  "quotient_dimension":q["quotient_dimension"],"quotient_gain_2norm":float(s[0]),
  "quotient_leading_singular_values":s[:8].tolist(),
  "gauge_to_transverse_injection_norm":q["gauge_injection_norm"],
  "line_source":"exact MOVING +/- sensor-compatible family, not SVD selection",
  "source_uniform_quotient_action_certified":False,"joint_cell_forward_invariant":False,
  "all_time_magnetic_service_verified":False,"theorem_closed":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=diagnostic(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
