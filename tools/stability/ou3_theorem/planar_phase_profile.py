"""Read and analyze OU3PRF1 native planar covariance profiles."""
from __future__ import annotations
import struct,json,numpy as np
from pathlib import Path

def read_profile(path):
 b=Path(path).read_bytes()
 if b[:7]!=b"OU3PRF1":raise ValueError("bad profile magic")
 n=struct.unpack_from("<I",b,8)[0];off=12;rows=[]
 stride=(2+144+81)
 vals=np.frombuffer(b,dtype="<f4",offset=off,count=n*stride).reshape(n,stride)
 for v in vals:
  rows.append((float(v[0]),float(v[1]),v[2:146].reshape(12,12).astype(float),v[146:].reshape(9,9).astype(float)))
 return rows

def candidate_radii(rows,inflation=1.25):
 if len(rows)!=8000:raise ValueError("need two 20-s periods")
 out={"even":[],"odd":[]}
 for k in range(4000):
  for name,i in (("even",2),("odd",3)):
   d=rows[k+4000][i]-rows[k][i]
   out[name].append(inflation*float(np.linalg.norm((d+d.T)/2,2)))
 return out

def summary(path):
 r=read_profile(path);rad=candidate_radii(r)
 return {"qualification":"OU3_PLANAR_PHASE_PROFILE_V1","samples":len(r),
         "period1_elapsed_start":r[0][0],"period2_elapsed_start":r[4000][0],
         "candidate_radius_even_max":max(rad["even"]),"candidate_radius_odd_max":max(rad["odd"]),
         "candidate_radius_inflation":1.25,"replay_radius_is_not_certificate":True,
         "containment_verified":False,"theorem_closed":False}
if __name__=="__main__":
 import sys;print(json.dumps(summary(sys.argv[1]),indent=2,sort_keys=True))
