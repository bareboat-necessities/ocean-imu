#!/usr/bin/env python3
"""Source-dynamic reachability backend for the OU-III P4 path certificate.

The old P4 worst cell combines low applied sigma_aw with very high applied R_S.
That corner is not fictitious: after a large-to-small sea transition the two
shipping EMA channels can lag by different amounts.  It is nevertheless a
*dynamical* source state, not an arbitrary Cartesian source choice.

This producer builds a conservative finite transition graph for the deployed
adaptation states (tau_applied, sigma_applied, R_S_applied) at the 0.1 s commit
cadence. tau/sigma share the deployed period-scaled EMA and R_S uses its own
1.5*tau target horizon. R_S_target is evaluated from the deployed default
SpectralMSE law using the same tau_target and sigma_target box. Edges are kept
whenever an outward-enlarged exact first-order EMA image intersects a destination
cell. No replay extrema are used.

The graph is a source-language certificate, not yet P4. It tells the next
complete-word Phi/Omega backend which weak cells can follow which other cells,
and whether the pathological old worst corner can persist indefinitely without
leaving that corner. This is exactly the correlation that the old global
min(delta) proof discarded.
"""
from __future__ import annotations

import argparse, json, math, re
from pathlib import Path

import ou3_source_reachable_matrix_p3 as P3
import ou3_source_domain_contract as SOURCE

REPO=Path(__file__).resolve().parents[1]
WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'

def down(x): return math.nextafter(float(x),-math.inf)
def up(x): return math.nextafter(float(x),math.inf)

def _literal_member(text,name):
    m=re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f",text)
    if not m: raise RuntimeError(f'cannot extract deployed literal member {name}')
    return float(m.group(1))

def _const(text,name): return float(SOURCE.parse_const(text,name))

def _constants():
    t=WRAPPER.read_text(encoding='utf-8')
    if 'RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;' not in t:
        raise RuntimeError('path backend requires deployed SpectralMSE law')
    return {
      'dt':_const(t,'FREQ_SMOOTHER_DT'),'commit':_const(t,'ADAPT_EVERY_SECS'),
      'tau_coeff':_literal_member(t,'tau_coeff_'),'sigma_coeff':_literal_member(t,'sigma_coeff_'),
      'adapt_tau_sea_periods':_const(t,'ADAPT_TAU_SEA_PERIODS'),'adapt_RS_mult':_const(t,'ADAPT_RS_MULT'),
      'min_tau':_const(t,'MIN_TAU_S'),'max_tau':_const(t,'MAX_TAU_S'),'max_sigma':_const(t,'MAX_SIGMA_A'),
      'min_RS':_const(t,'MIN_R_S'),'max_RS':_const(t,'MAX_R_S'),
      'min_freq':_const(t,'MIN_TUNE_FREQ_HZ'),'max_freq':_const(t,'MAX_TUNE_FREQ_HZ'),
      'pseudo_ratio':_const(t,'PSEUDO_UPDATE_TAU_RATIO_DEFAULT'),
      'pseudo_min':_const(t,'PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT'),'pseudo_max':_const(t,'PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT'),
      'mse_coeff':_const(t,'R_S_MSE_COEFF_DEFAULT'),'noise_density':_const(t,'R_S_ACCEL_NOISE_DENSITY_DEFAULT'),
    }

def _target_rs(tau,sigma,c):
    tau=min(max(float(tau),c['min_tau']),c['max_tau']); sigma=min(max(float(sigma),1e-6),c['max_sigma'])
    TS=min(max(c['pseudo_ratio']*tau,c['pseudo_min']),c['pseudo_max'])
    qpow=(2*c['noise_density'])**(1/14); sigma_aB=sigma/c['sigma_coeff']
    raw=c['mse_coeff']*qpow*(sigma_aB**(6/7))*(tau**(24/7))/math.sqrt(TS)
    return min(max(raw,c['min_RS']),c['max_RS'])

def _rs_box(tau,sigma,c): return (down(_target_rs(tau[0],sigma[0],c)),up(_target_rs(tau[1],sigma[1],c)))

def _ema_image(x,target,horizon,dt):
    if not horizon[0]>0 or horizon[1]<horizon[0]: raise RuntimeError('invalid EMA horizon')
    alo=down(math.exp(-dt/horizon[0])); ahi=up(math.exp(-dt/horizon[1])); vals=[]
    for a in (alo,ahi):
      for xx in x:
       for u in target: vals.append(a*xx+(1-a)*u)
    return (down(min(vals)),up(max(vals)))

def _overlap(a,b): return not (a[1]<b[0] or b[1]<a[0])
def _cells(lo,hi,n):
    e=P3.geom_edges(float(lo),float(hi),int(n)); return [(down(e[i]),up(e[i+1])) for i in range(len(e)-1)]
def _matching(cells,image): return [i for i,c in enumerate(cells) if _overlap(c,image)]
def _tau_target(freq,c):
    lo=c['tau_coeff']*.5/freq[1]; hi=c['tau_coeff']*.5/freq[0]
    return (down(max(c['min_tau'],min(c['max_tau'],lo))),up(max(c['min_tau'],min(c['max_tau'],hi))))

def _scc(graph):
    n=len(graph); seen=[False]*n; order=[]
    def dfs(v):
      seen[v]=True
      for w in graph[v]:
       if not seen[w]: dfs(w)
      order.append(v)
    for v in range(n):
      if not seen[v]: dfs(v)
    rg=[[] for _ in range(n)]
    for v,ws in enumerate(graph):
      for w in ws: rg[w].append(v)
    seen=[False]*n; comps=[]
    def rdfs(v,a):
      seen[v]=True;a.append(v)
      for w in rg[v]:
       if not seen[w]: rdfs(w,a)
    for v in reversed(order):
      if not seen[v]: a=[];rdfs(v,a);comps.append(a)
    return comps

def _induced_cycle(graph,nodes):
    ns=set(nodes); sub={v:[w for w in graph[v] if w in ns] for v in ns}; visiting=set();done=set()
    def dfs(v):
      if v in visiting:return True
      if v in done:return False
      visiting.add(v)
      for w in sub[v]:
       if dfs(w):return True
      visiting.remove(v);done.add(v);return False
    return any(dfs(v) for v in list(ns) if v not in done)

def _longest_bad_residence(graph,bad):
    ns=set(bad); sub={v:[w for w in graph[v] if w in ns] for v in ns}
    if _induced_cycle(graph,bad): return None
    memo={}
    def f(v):
      if v not in memo:memo[v]=1+max((f(w) for w in sub[v]),default=0)
      return memo[v]
    return max((f(v) for v in ns),default=0)

def build(domain_path=DEFAULT_DOMAIN):
    dom=json.loads(Path(domain_path).read_text(encoding='utf-8'))
    if dom.get('trajectory_fit') is not False: raise RuntimeError('path domain must not be trajectory fitted')
    c=_constants(); tau=_cells(max(c['min_tau'],c['tau_coeff']*.5/c['max_freq']),c['max_tau'],10); sig=_cells(.05,c['max_sigma'],8); rs=_cells(c['min_RS'],c['max_RS'],10); freq=_cells(c['min_freq'],c['max_freq'],8)
    states=[]; idx={}
    for i,t in enumerate(tau):
      for j,s in enumerate(sig):
       for k,r in enumerate(rs): idx[(i,j,k)]=len(states);states.append((t,s,r))
    targets=[]
    for f in freq:
      tt=_tau_target(f,c)
      for ss in sig:
       rr=_rs_box(tt,ss,c); ht=(down(c['adapt_tau_sea_periods']*tt[0]/c['tau_coeff']),up(c['adapt_tau_sea_periods']*tt[1]/c['tau_coeff'])); hr=(down(c['adapt_RS_mult']*tt[0]),up(c['adapt_RS_mult']*tt[1])); targets.append((tt,ss,rr,ht,hr))
    graph=[set() for _ in states]; dt=c['commit']
    for q,(t,s,r) in enumerate(states):
      out=graph[q]
      for tt,ss,rr,ht,hr in targets:
       ti=_matching(tau,_ema_image(t,tt,ht,dt)); si=_matching(sig,_ema_image(s,ss,ht,dt)); ri=_matching(rs,_ema_image(r,rr,hr,dt))
       for i in ti:
        for j in si:
         for k in ri: out.add(idx[(i,j,k)])
    gl=[sorted(x) for x in graph]; comps=_scc(gl); recurrent=set()
    for cc in comps:
      if len(cc)>1 or (cc and cc[0] in graph[cc[0]]): recurrent.update(cc)
    bad=[]
    for q,(t,s,r) in enumerate(states):
      x=(c['dt']/t[1],c['dt']/t[0])
      if _overlap(s,(.05,.13025855423486765)) and _overlap(r,(149.21548743644342,400.0)) and _overlap(x,(.00041666665735344083,.0004837652693428343)): bad.append(q)
    bad_cycle=_induced_cycle(gl,bad); steps=_longest_bad_residence(gl,bad)
    return {
      'qualification':'OU3_P4_SOURCE_DYNAMIC_PATH_REACHABILITY','source_only':True,'trajectory_replay_used':False,'deployed_default_law':'SpectralMSE','commit_period_s':dt,
      'partition':{'tau':len(tau),'sigma':len(sig),'R_S':len(rs),'states':len(states),'target_boxes':len(targets)},'transition_edges':sum(map(len,gl)),'strongly_connected_components':len(comps),'recurrent_states':len(recurrent),
      'old_worst_corner_state_count':len(bad),'old_worst_corner_states_in_any_recurrent_SCC':sum(q in recurrent for q in bad),'old_worst_corner_has_internal_recurrent_cycle':bad_cycle,
      'old_worst_corner_max_consecutive_commit_steps_upper':steps,'old_worst_corner_max_residence_s_upper':None if steps is None else up(steps*dt),
      'path_graph_ready':True,'usable_P4_promoted':False,'next_obligation':'propagate complete-word Phi/Omega and exact nonlinear return map on this source-reachable graph; charge weak cells by reachable residence/path products instead of global min(delta)','failures':[]}

def validate(d):
    f=list(d.get('failures',[]))
    if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False:f.append('path graph is not source-only')
    if d.get('path_graph_ready') is not True:f.append('path graph not ready')
    if int(d.get('partition',{}).get('states',0))<=0 or int(d.get('transition_edges',0))<=0:f.append('empty path graph')
    if d.get('usable_P4_promoted') is not False:f.append('reachability prematurely promoted P4')
    if int(d.get('old_worst_corner_state_count',0))<=0:f.append('old worst corner not represented')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain.resolve());f=validate(d);d['validation_failures']=f;a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding='utf-8');print(json.dumps({k:d[k] for k in ('partition','transition_edges','strongly_connected_components','recurrent_states','old_worst_corner_state_count','old_worst_corner_states_in_any_recurrent_SCC','old_worst_corner_has_internal_recurrent_cycle','old_worst_corner_max_residence_s_upper','failures')},indent=2,sort_keys=True));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())
