#!/usr/bin/env python3
"""Reconstruct the P3 limiting translation source corner without a 2000-cell scan.

The current P4 diagnostic showed that the old source-uniform margin is limited
by the corner with smallest h/tau, smallest sigma_aw, and largest R_S.  This
helper reconstructs that exact partition cell from the same source-derived P3
grids and evaluates only that cell with the direct generalized-matrix backend.
It is used for rapid iteration of the complete-word replacement proof.

This is not a new assumption: `verify_against_global_scan()` can be used in CI
or audit code to compare the reconstructed corner with the exhaustive P3 scan.
"""
from __future__ import annotations

import json, math
from pathlib import Path

import ou3_source_reachable_matrix_p3_direct as P3D

REPO=Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'


def build_cell(mode: str, domain_path: Path=DEFAULT_DOMAIN) -> dict:
 p=Path(domain_path).resolve(); domain=json.loads(p.read_text(encoding='utf-8'))
 live=domain['normal_live']; B=P3D.BASE
 vector=B.VECTOR.build(); process=B.PROCESS.build(); words=B.WORDS.build(p)
 for label,obj,val in (('vector',vector,B.VECTOR.validate),('process',process,B.PROCESS.validate),('words',words,B.WORDS.validate)):
  f=val(obj)
  if f: raise RuntimeError(f'{label} prerequisite failed: {f}')
 alpha6=B.vector_alpha6(live,vector); sched=B.source_schedule(); h=float(sched['dt_s'])
 tau_lo,tau_hi=map(float,sched['tau_applied_invariant_s']); xlo,xhi=h/tau_hi,h/tau_lo
 edges=B.geom_edges(xlo,xhi,24)
 if xlo<B.BRANCH_X<xhi: edges=sorted(set(edges+[B.BRANCH_X]))
 xcells=[]
 for c in B.interval_cells(edges): xcells.extend(B.split_x_cell(c))
 sigmas=B.interval_cells(B.geom_edges(0.05,6.0,5))
 rs_lo,rs_hi=sched['R_S_applied_invariant']; rss=B.interval_cells(B.geom_edges(rs_lo,rs_hi,8))
 x,rho=xcells[0]; sigma=sigmas[0]; rs=rss[-1]
 row=P3D.mode_cell(mode,x,rho,sigma,rs,live,vector,process,sched,alpha6)
 g=row['generalized_matrix_inequality']
 if g['limiting_block']!='translation_RL_inverse_block':
  raise RuntimeError(f'{mode}: reconstructed limiting corner changed: {g["limiting_block"]}')
 return {'mode':mode,'x':x,'rho_translation_lower':rho,'sigma':sigma,'rs':rs,'row':row,'sched':sched,'vector':vector,'process':process,'live':live}


def serializable(c:dict)->dict:
 r=c['row']
 return {'mode':c['mode'],'x_h_over_tau':c['x'].as_list(),'sigma_aw':c['sigma'].as_list(),'R_S':c['rs'].as_list(),'rho_translation_lower':float(c['rho_translation_lower']),'delta_full_lower':float(r['relative_Riccati_injection_margin_lower']),'delta_translation_lower':float(r['direct_translation_generalized_margin_lower']),'delta_nontranslation_lower':float(r['direct_nontranslation_margin_lower']),'limiting_block':r['generalized_matrix_inequality']['limiting_block']}
