"""80-digit linked-supply diagnostic on carried source words; no promotion.

Reports prospective local-defect parity and endpoint source attribution before any enclosure.
Shipping-faithfulness protocol applies: these carried diagnostics are non-promoting.

Use --export to compile an observer and untapped control via the existing
readout fixture. Otherwise analyze previously exported word-{0,wave}.json.
Actual endpoint defects are retrospective; they are not uniform supply bounds.
"""
from pathlib import Path
import argparse
import json
import hashlib
import shutil
import subprocess
import mpmath as mp

from .ag_readout_source_diagnostic import instrument, HEADER
REPO = Path(__file__).resolve().parents[3]

def mat(a): return mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
def sym(a): return (a+a.T)/2
def skew(v):
 x,y,z=v
 return mp.matrix([[0,-z,y],[z,0,-x],[-y,x,0]])
def qprod(a,b):
 va,vb=mp.matrix(a[:3]),mp.matrix(b[:3]); v=a[3]*vb+b[3]*va+skew(va)*vb
 return list(v)+[a[3]*b[3]-(va.T*vb)[0]]
def error(state, qraw, t, live, wave):
 # NED, error q_true=Exp(delta_theta)q_hat (left W->B error).
 q=[mp.mpf(float(row[0])) for row in qraw]; norm=mp.sqrt(sum(x*x for x in q));q=[x/norm for x in q]
 angle=-mp.mpf('.02')*mp.sin(t/2) if wave else mp.mpf(0)
 qt=[mp.sin(angle/2),mp.mpf(0),mp.mpf(0),mp.cos(angle/2)]
 qe=qprod(qt,[-q[0],-q[1],-q[2],q[3]])
 if qe[3]<0: qe=[-x for x in qe]
 n=mp.sqrt(sum(x*x for x in qe[:3]));rv=mp.matrix(qe[:3])*(2*mp.atan2(n,qe[3])/n) if n else mp.zeros(3,1)
 e=-mat(state);e[:3,0]=rv
 if wave:
  e[8]+=mp.mpf('.24')*mp.cos(mp.mpf('.6')*t)
  e[11]+=mp.mpf('.4')*mp.sin(mp.mpf('.6')*t)
  e[14]+=mp.mpf(2)/3*(mp.cos(mp.mpf('.6')*live)-mp.cos(mp.mpf('.6')*t))
  e[17]-=mp.mpf('.144')*mp.sin(mp.mpf('.6')*t)
 return e
fmt=lambda x:mp.nstr(x,25)
def analyze(path):
 trace=json.loads(path.read_text())
 M=mp.eye(21)
 for ev in trace['events']:
  if ev['kind']=='prediction':
   M[:6,:]=mat(ev['F_AG'])*M[:6,:]
   M[6:18,:]=mat(ev['F_LIN'])*M[6:18,:]
   M[18:,:]=mp.mpf(float(ev['phi_BA']))*M[18:,:]
  elif ev['kind']=='correction': M=M-mat(ev['K'])*(mat(ev['H'])*M)
  elif ev['kind']=='reset':
   v=[mp.mpf(float(x[0])) for x in ev['d']]
   M[:3,:]=(mp.eye(3)+skew(v)/2)*M[:3,:]
  elif ev['kind'] in ('sync','sync_completion','adaptive_state'): pass
  else: raise ValueError(ev['kind'])
 P0=sym(mat(trace['root_covariance'])); PN=sym(mat(trace['terminal_covariance']))
 L0=mp.cholesky(P0);LN=mp.cholesky(PN); Li=LN**-1
 H=Li*M*L0
 Delta=sym(mp.eye(21)-H.T*H); eig=mp.eigsy(Delta,eigvals_only=True); dmin=min(eig)
 from .local_defect_composition import carried_boundaries
 B=carried_boundaries(trace['events'])
 e0=mat(B[0]['e_before']); e1=mat(B[-1]['e_after'])
 # Prospective 80-digit composition for long carried diagnostics.  The exact
 # Fraction implementation remains the algebraic certificate; using it for a
 # 100-s/20k-step feasibility replay creates enormous rational numerators with
 # no theorem benefit.
 b_local=mp.zeros(21,1); M_local=mp.eye(21); grouped={}
 for item in B:
  A=mat([[float(v) for v in row] for row in item['A']])
  eb=mat(item['e_before']); ea=mat(item['e_after'])
  d=ea-A*eb
  for key in list(grouped): grouped[key]=A*grouped[key]
  label=item['kind'] if item['kind']!='correction' else 'correction:'+str(item.get('sensor'))
  grouped[label]=grouped.get(label,mp.zeros(21,1))+d
  b_local=A*b_local+d
  M_local=A*M_local
 prospective_parity=mp.norm(e1-M_local*e0-b_local)
 homogeneous_map_parity=mp.norm(M-M_local)
 x=L0**-1*e0; b=Li*b_local; z=H.T*b
 grouped_metric={key:fmt(mp.norm(Li*value)) for key,value in grouped.items()}
 block_names=(('theta',0,3),('bg',3,6),('v',6,9),('p',9,12),('S',12,15),('aw',15,18),('ba',18,21))
 endpoint_raw_blocks={name:fmt(mp.norm(b_local[a:z0])) for name,a,z0 in block_names}
 V0=(x.T*x)[0]; VN=(e1.T*(PN**-1)*e1)[0]
 ans={'profile':trace['profile'],'dps':80,'event_count':len(trace['events']),
  'observer_terminal_parity':trace['observer_terminal_parity'],
  'trace_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
  'source_header_sha256':trace['source_header_sha256'],
  'root_sqrt_V':fmt(mp.sqrt(V0)),'terminal_sqrt_V':fmt(mp.sqrt(VN)),
  'V_change':fmt(VN-V0),'homogeneous_rho':fmt(1-dmin),
  'minimum_whitened_loss':fmt(dmin),'signed_endpoint_forcing_norm':fmt(mp.norm(b)),
  'active_release_time_s':trace['active_step']/200,
  'live_time_s':trace['live_step']/200,
  'release_LIN_mean_norm':fmt(mp.norm(mat(trace['release_state'])[6:18,:])),
  'source_uniform_verified':False,'theorem_closed':False,
  'supply_is_retrospective_endpoint_residual':False,
  'local_defect_origin':'prospective_literal_same_boundary_composition',
  'local_b_endpoint_parity_norm':fmt(prospective_parity),
  'homogeneous_map_parity_norm':fmt(homogeneous_map_parity),
  'endpoint_forcing_metric_norm_by_operation_class':grouped_metric,
  'structures_preserved':['literal mean chronology','literal covariance-derived gains','persistent physical history','persistent filter mean/covariance','shipping S scheduler and covariance sync chronology'],
  'relaxations_introduced':['finite carried fixture is diagnostic only; not source-uniform','scalar chi/gamma remains a sufficient-bound diagnostic, not a shipping impossibility test'],
  'failure_class_if_ratio_misses_target':'D_SUFFICIENT_BOUND_FAILURE',
  'endpoint_forcing_raw_norm_by_state_block':endpoint_raw_blocks,
  'uniform_inner_retention_refuted_by_this_finite_replay':False,
  'all_time_magnetic_service_verified':False,
  'physical_S_origin':'one fixed capture/Live epoch; never reset at word boundaries',
  'metric':'exact symmetric part of exported float covariance; arithmetic enclosure not claimed'}
 if dmin<=0: ans['status']='NO_POSITIVE_FROZEN_WORD_LOSS';return ans
 def eval_gamma(g):
  GG=Delta-g*mp.eye(21)
  cc=(b.T*b)[0]+(z.T*(GG**-1)*z)[0]
  return cc/g,cc,GG
 # Exact scalar optimality condition for f(gamma)=chi(gamma)/gamma:
 # gamma*chi'(gamma)-chi(gamma)=0.  Since chi''>=0, its left side is
 # monotone nondecreasing; bisection therefore finds the unique interior
 # minimizer when it exists.  A strict-boundary infimum is approached from
 # below and is sufficient for the diagnostic.
 c0=(b.T*b)[0]
 def stationarity(g):
  GG=Delta-g*mp.eye(21)
  inv=GG**-1
  chi=c0+(z.T*inv*z)[0]
  chip=(z.T*inv*inv*z)[0]
  return g*chip-chi
 eps=dmin*mp.mpf('1e-30')
 lo=dmin*mp.mpf('1e-30'); hi=dmin-eps
 if c0==0 and mp.norm(z)==0:
  gam=dmin/2
  optimizer_case='zero_forcing'
 elif stationarity(hi)<=0:
  gam=hi
  optimizer_case='boundary_infimum'
 else:
  a,c=lo,hi
  for _ in range(240):
   mid=(a+c)/2
   if stationarity(mid)>0: c=mid
   else: a=mid
  gam=(a+c)/2
  optimizer_case='unique_interior_stationary'
 ratio,chi,G=eval_gamma(gam)
 completed=(1-gam)*V0+chi-((x-G**-1*z).T*G*(x-G**-1*z))[0]
 ans.update({'defect_origin':'prospective_local_defect_composition',
  'eligible_for_theorem_entry_test':False,
  'gamma':fmt(gam),'gamma_optimizer_case':optimizer_case,'linked_supply_chi':fmt(chi),
  'optimized_chi_over_gamma':fmt(ratio),
  'optimized_linked_to_inner_budget_ratio':fmt(ratio/(mp.mpf('.15')**2)),
  'frozen_fixed_forcing_sufficient_radius':fmt(mp.sqrt(ratio)),
  'linked_to_inner_budget_ratio':fmt(chi/(gam*mp.mpf('.15')**2)),
  'linked_to_root_budget_ratio':fmt(chi/(gam*V0)) if V0 else None,
  'separated_fixed_forcing_radius':fmt(mp.norm(b)/(1-mp.sqrt(1-dmin))),
  'completed_square_identity_error':fmt(abs(VN-completed)),
  'status':'NON_PROMOTING_LINKED_FEASIBILITY_ONLY'})
 return ans


def export_words(directory, eigen, compiler, horizons=(16,30,60,100)):
    directory.mkdir(parents=True, exist_ok=True)
    inc = directory/'include'/'kalman_ou_iii'
    inc.mkdir(parents=True, exist_ok=True)
    header = (REPO/HEADER).read_text()
    (inc/HEADER.name).write_text(instrument(header))
    source0 = (REPO/'tools/stability/ag_readout_source.cpp').read_text()
    for horizon in horizons:
        source = source0
        def once(old, new):
            nonlocal source
            if source.count(old) != 1:
                raise ValueError('readout driver anchor changed: '+old[:70])
            source = source.replace(old, new)
        samples=int(round(float(horizon)/.005))
        terminal=45000+samples
        once('for (int k=1; k<=45064; ++k)',f'for (int k=1; k<={terminal}; ++k)')
        once('    std::string root;', '    std::string root, root_state, root_quat, release_state;')
        once('            recording=true;',
             '            root_state=matrix_json(filter.raw().mekf().xext);\n'
             '            root_quat=matrix_json(filter.raw().mekf().qref.coeffs());\n'
             '            recording=true;')
        once('if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;',
             'if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) {\n'
             '          active=k; release_state=matrix_json(filter.raw().mekf().xext);\n        }')
        once('if (!m.Pext.allFinite() || !m.xext.allFinite() || applied!=8) return 3;',
             'if (!m.Pext.allFinite() || !m.xext.allFinite()) return 3;')
        needle = '              << ",\\\"terminal_covariance\\\":" << matrix_json(m.Pext)'
        once(needle,
             '              << ",\\\"root_state\\\":" << root_state\n'
             '              << ",\\\"root_quaternion\\\":" << root_quat\n'
             '              << ",\\\"release_state\\\":" << release_state\n'+needle)
        driver = directory/f'linked-source-{horizon}s.cpp'
        driver.write_text(source)
        source_hashes = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted((REPO/'src').rglob('*.h'))}
        binaries={}
        for name, includes in [('observed', ['-I'+str(inc.parent)]), ('control', [])]:
            binary=directory/f'{name}-{horizon}s'
            subprocess.run([compiler, '-O1', '-std=c++20', *includes,
                            '-I'+str(REPO/'src'), '-isystem', str(eigen),
                            str(driver), '-o', str(binary)], check=True)
            binaries[name]=binary
        for profile in ['0', 'wave']:
            obs = json.loads(subprocess.check_output([str(binaries['observed']), profile], text=True))
            ctl = json.loads(subprocess.check_output([str(binaries['control']), profile], text=True))
            if any(obs[k] != v for k, v in ctl.items() if k != 'events'):
                raise ArithmeticError('observer changed source execution')
            obs.update(profile=profile, horizon_s=horizon, observer_terminal_parity=True,
                       source_header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                       source_hashes=source_hashes)
            (directory/f'word-{profile}-{horizon}s.json').write_text(json.dumps(obs, sort_keys=True))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--export', action='store_true')
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--cxx', default=shutil.which('clang++') or 'g++')
    parser.add_argument('--horizons', default='16,30,60,100')
    args = parser.parse_args()
    horizons=tuple(int(x) for x in args.horizons.split(','))
    if args.export:
        export_words(args.directory.resolve(), args.eigen, args.cxx, horizons)
    with mp.workdps(80):
        report = {'kind': 'NON_PROMOTING_LINKED_FINITE_SUPPLY', 'decimal_digits': 80,
                  'source_uniform_verified': False, 'theorem_closed': False,
                  'words': [analyze(args.directory/f'word-{p}-{h}s.json')
                            for h in horizons for p in ['0','wave']]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\\n')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
