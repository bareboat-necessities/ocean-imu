"""Generate C++ initializer for certified anisotropic 21x21 root lower matrix."""
from __future__ import annotations
from .planar_anisotropic_factors import full_root_lower

def cpp_initializer():
 p=full_root_lower()
 vals=[]
 for j in range(21):
  for i in range(21): vals.append(f"{p[i,j]:.17g}")
 return "{" + ",".join(vals) + "}"

if __name__=="__main__": print(cpp_initializer())
