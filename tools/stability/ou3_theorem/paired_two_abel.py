"""Paired two-Abel physical-primitive to causal-source certificate layer."""
from fractions import Fraction as F
def collapse_slab_boundaries(endpoint_rows):
 rows=[([F(x) for x in l],[F(x) for x in r]) for l,r in endpoint_rows]
 if not rows: raise ValueError("slabs required")
 out=[rows[0][0]]
 for (_,r),(l,_) in zip(rows[:-1],rows[1:]):
  if len(r)!=len(l): raise ValueError("matching dimensions")
  out.append([a+b for a,b in zip(r,l)])
 out.append(rows[-1][1]); return out
def paired_source_columns(process_columns,accel_columns):
 if len(process_columns)!=len(accel_columns): raise ValueError("matching rows")
 out=[]
 for p,a in zip(process_columns,accel_columns):
  if len(p)!=len(a): raise ValueError("matching columns")
  out.append([F(x)+F(y) for x,y in zip(p,a)])
 return out
def frobenius_sq(m): return sum((F(x)*F(x) for r in m for x in r),F(0))
def certificate_schema():
 return {"qualification":"OU3_PAIRED_TWO_ABEL_V1","uses_C_port_W":False,
 "pairs_process_and_accelerometer_before_norm":True,
 "telescopes_internal_physical_boundaries_before_norm":True,
 "requires_literal_F_Q_K_H_S_sync_reset_order":True,
 "requires_same_history_physical_p_v_S_a_lift":True,
 "fixed_word_exact_algebra_certificate_available":True,
 "outer_class_factor_interval_enclosure_available":False,
 "source_uniform_induced_norm_numeric":None,"candidate_reserve_rad":.003056118,
 "strict_bridge_closed":False,"theorem_closed":False}
