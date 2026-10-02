"""Cellwise rigorous s_i/d_i composition for reachable shaped storage."""
from __future__ import annotations
import math
from dataclasses import dataclass
from .interval_riccati import symmetric_interval_gershgorin

@dataclass(frozen=True)
class CellRatio:
 token:str;d_lower:float;s_upper:float;ratio_upper:float;verified:bool;reason:str

def homogeneous_lower(iq,kernel_certificate:dict):
 if not kernel_certificate.get("zero_set_excluded",False):
  return 0.,False,"zero set not excluded"
 # This lower bound is on the already kernel-restricted history cell. The
 # caller must provide a restriction certificate; full-space Gershgorin alone
 # is not allowed to claim kernel removal.
 R=kernel_certificate.get("restricted_action")
 if R is None:return 0.,False,"kernel-restricted action missing"
 lo,hi=symmetric_interval_gershgorin(R.mid,R.rad)
 return max(0.,lo),lo>0,"ok" if lo>0 else "restricted action not positive"

def combine(token,iq,kernel_certificate,supply_certificate):
 d,ok,reason=homogeneous_lower(iq,kernel_certificate)
 if not ok:return CellRatio(token,d,math.inf,math.inf,False,reason)
 if not supply_certificate.verified:return CellRatio(token,d,supply_certificate.upper,math.inf,False,supply_certificate.reason)
 s=max(0.,supply_certificate.upper);return CellRatio(token,d,s,s/d,True,"verified")

def cover_max(cells):
 if not cells:return {"verified":False,"reason":"empty cover"}
 bad=[x for x in cells if not x.verified]
 if bad:return {"verified":False,"reason":"unverified cell","first_bad":bad[0]}
 m=max(cells,key=lambda x:x.ratio_upper)
 return {"verified":True,"max_ratio_upper":m.ratio_upper,"limiting_token":m.token,
         "cell_count":len(cells)}

def kernel_certificate_from_words(A0,line0,T,A1,line1):
 from .kernel_restricted_action import certify_rotating_line
 R,c=certify_rotating_line(A0,line0,T,A1,line1)
 return {"zero_set_excluded":bool(c["verified"]),"restricted_action":R,
         "restricted_eigen_lower":c["restricted_eigen_lower"],
         "word_dependent_kernel":True,"later_word_action_included":True,
         "line_rotation_rad":c["midpoint_line_rotation_rad"]}

def kernel_certificate_from_line_cone(A0,line0,line0_radius,T,A1,line1):
 from .kernel_restricted_action import augment_later_word_action,restrict_interval_action
 As=augment_later_word_action(A0,T,A1)
 R,c=restrict_interval_action(As,line0,line0_radius)
 return {"zero_set_excluded":bool(c["verified"]),"restricted_action":R,
         "restricted_eigen_lower":c["restricted_eigen_lower"],
         "word_dependent_kernel":True,"later_word_action_included":True,
         "line_cone_covered":c.get("whole_line_cone_covered",False),
         "line_angle_rad":c.get("line_angle_rad",0.)}
