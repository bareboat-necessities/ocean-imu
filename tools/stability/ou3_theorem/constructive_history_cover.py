"""Executable constructive history-cover pipeline.

Consumes causal leaf objects produced by the shipping history propagator.
Promotion is all-or-nothing: every leaf must verify nonlinear premises,
kernel-restricted action, and direct QCQP supply.
"""
from __future__ import annotations
from dataclasses import dataclass
import json,math
from pathlib import Path
import numpy as np
from .verified_linked_qcqp import Box,verified_interval_supply_bnb,shipping_side_callback
from .constructive_cell_ratio import kernel_certificate_from_line_cone,combine,cover_max
from .marine_magnetic_qcqp import callback
ROOT=Path(__file__).resolve().parents[3]
CONSTANTS=json.loads((ROOT/"tools/stability/ou3_theorem/constants.json").read_text())

@dataclass
class ConstructiveLeaf:
 token:str;witness:object;interval_joint:object
 A0:object;T:object;A1:object
 e_box:Box;u_box:Box;linear_A:np.ndarray;linear_b:np.ndarray

def certify_leaf(x:ConstructiveLeaf,*,tol=1e-4,max_leaves=200000):
 if callback(x.witness.marine_magnetic,CONSTANTS) is not True:
  return {"verified":False,"token":x.token,"reason":"nonlinear MARINE/MAGNETIC unresolved"}
 kc=kernel_certificate_from_line_cone(x.A0,x.witness.kernel_line.midpoint,
       x.witness.kernel_line.component_radius,x.T,x.A1,x.witness.kernel_line.midpoint)
 side=shipping_side_callback(lambda box:x.witness.marine_magnetic,CONSTANTS)
 supply=verified_interval_supply_bnb(x.interval_joint.B,x.interval_joint.C,
       x.e_box,x.u_box,x.linear_A,x.linear_b,side,tol=tol,max_leaves=max_leaves)
 z=combine(x.token,x.interval_joint,kc,supply)
 return {"verified":z.verified,"token":x.token,"ratio":z,"kernel":kc,"supply":supply}

def certify_cover(leaves):
 out=[certify_leaf(x) for x in leaves]
 if any(not x["verified"] for x in out):
  return {"verified":False,"leaves":out,"reason":"first unverified leaf"}
 ratios=[x["ratio"] for x in out]
 return {"verified":True,"max":cover_max(ratios),"leaves":out}
